
/*!
 /** Copyright (c) CDC IXIS CM July 2003 Paris
 *	\file CRAcalculator.cpp
 *  \brief file for the Bermuda Swaption Calculator
 *	\author  H. BAKHTRI & M. ABDELMOUMNI & P. LAM
 *	\version 1.0
 *	\date June 2005
 */

/// this header comes first as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpcalculators/cracalculator.h"
#include "gpcalculators/argconvdefault.h"

/// gpbase
#include "gpbase/env.h"
#include "gpbase/typedef.h"
#include "gpbase/argconvdefault.h"
#include "gpbase/datestripcombiner.h"
#include "gpbase/autocleaner.h"
#include "gpbase/singleton.h"
#include "gpbase/datestrip.h"
#include "gpbase/stringconvert.h"
#include "gpbase/gptrigomatrix.h"
#include "gpbase/gpvector.h"
#include "gpbase/curveconvert.h"

#include "gpbase/ostringstream.h"
#include "gpbase/curve.h"
#include "gpbase/gplinalgconvert.h"
#include "gpbase/numericconstant.h"
#include "gpbase/surface.h"
#include "gpbase/surfacetypedef.h"

/// gpinfra
#include "gpinfra/argconvdefault.h"
#include "gpinfra/gensecurity.h"
#include "gpinfra/dealdescription.h"
#include "gpinfra/gensecmanipulator.h"
#include "gpinfra/mktdatamanagerrep.h"
#include "gpinfra/pricingadviser.h"
#include "gpinfra/genpricer.h"
#include "gpinfra/numeraire.h"
#include "gpinfra/numerairefactory.h"
#include "gpinfra/pricerinfo.h"
#include "gpinfra/gramfunctorargdict.h"
#include "gpinfra/cstmanager.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/correlmatparam.h"
#include "gpinfra/modelnamemap.h"

/// gpcalib
#include "gpcalib/modelparamsfactory.h"
#include "gpcalib/calibmethod.h"
#include "gpcalib/vanillacap.h"
#include "gpcalib/kerneltogp.h"
#include "gpcalib/modelfitter.h"
#include "gpcalib/modelfitterdes.h"
#include "gpcalib/vanillapricer.h"

/// gpmodels
#include "gpmodels/modelparamshw1f.h"
#include "gpmodels/sfrm.h"
#include "gpmodels/modelparamssfrm.h"
#include "gpmodels/modelparamssfrmfactory.h"
#include "gpmodels/MultiAssetsFactory.h"
#include "gpmodels/multiassets.h"
#include "gpmodels/local_SLN_Model.h"
#include "gpmodels/local_SLN_ModelParams.h"

/// gpnumlib
#include "gpnumlib/randomgenfactory.h"
#include "gpnumlib/compositegen.h"
#include "gpnumlib/antitheticgen.h"
#include "gpnumlib/argconvdefault.h"

/// gpnummethods
/// Tree
#include "gpnummethods/treemethod.h"
#include "gpnummethods/treebase.h"
#include "gpnummethods/treefactory.h"
#include "gpnummethods/scheduler.h"
#include "gpnummethods/sampler.h"
#include "gpnummethods/truncator.h"
#include "gpnummethods/reconnector.h"
#include "gpnummethods/smoother.h"
/// AMC
#include "gpnummethods/amcmethod.h"
#include "gpnummethods/amc_andersen.h"
#include "gpnummethods/amc_ls.h"
#include "gpnummethods/cfmethod.h"

/// kernel
#include "inst/swaption.h"
#include "inst/corridorleg.h"
#include "inst/swapleg.h"
#include "inst/fixleg.h"
#include "inst/portfolio.h"
#include "glob/paramview.h"
#include "zeroflat.h"
#include "volcube.h"
#include "volflat.h"
#include "armfrmhwvol.h"

/// STL
#include <iomanip> /// for setprecision()
#include <memory>

CC_BEGIN_NAMESPACE( ARM )

/// Default MDM key names
const string YC_KEY_NAME            = "YC_";
const string OSWMODEL_KEY_NAME      = "OSWMOD_";
const string CFMODEL_KEY_NAME       = "CFMOD_";
const string MRS_KEY_NAME           = "MRS_";
const string BETA_KEY_NAME          = "BETA_";

/// Call schedule 
const unsigned int CALL_CRA_SCHED	= 0;
const double NON_CALL_FEE			= 1.0e15;

/// SFRM sigma range [10bp,1000bp] with a 20% default value
const double SIGMA_LOWER_BOUND      = 0.001;
const double SIGMA_UPPER_BOUND      = 1.0;
const double SIGMA_DEFAULT_VALUE    = 0.01;

/// SFRM BETA range [1%,150%] with a 100% default value
const double BETA_LOWER_BOUND       = 0.01;
const double BETA_UPPER_BOUND       = 1.5;

const int SFRM_NB_FACTORS			= 1;

const double treeParams[]			= {3.0, 30.0, 7.0, 2.0, 0.0};
bool skipCalib						= false;

const int DEAL_NBCOLUMNS			= 18;
const int DEAL_FLOOR_NBCOLUMNS		= 17;
const int DEAL_CAPFLOOR_NBCOLUMNS	= 13;


const string ARM_CRACalculator::CRAColNamesTable [] =
{  
	"ResetDate",
	"StartDate",
	"NextStartDate",
	"MaturityDate",
	"Fees",
	"Funding",
	"Corridor",
	"ExoticSwap",
	"ExoticSwap1",
	"Option",  
	"Bermuda",
	"FundingLeg",
	"CorridorLeg",
	"ExoticSwap2",
	"Option2", 
	"Bermuda2",
	"FwdCorridor",
	"FwdFunding",
};


const int ARM_CRACalculator::CRAProductToPriceColumns [] =
{
    Bermuda,
	Corridor,
    Funding,
	FwdCorridor,
	FwdFunding,
    Bermuda2,
	ExoticSwap,
};

////////////////////////////////////////////////////
///	Class  : ARM_CRACalculator
///	Routine: CopyNoCleanUp
///	Returns: void
///	Action : Arguments copy
////////////////////////////////////////////////////
void ARM_CRACalculator::CopyNoCleanUp(const ARM_CRACalculator& rhs)
{
	itsCcy						= (rhs.itsCcy);
	itsStartDate				= (rhs.itsStartDate);
	itsEndDate					= (rhs.itsEndDate);
	itsPayRec					= (rhs.itsPayRec);
	itsNotional					= (rhs.itsNotional);
	itsRealFundNotional         = rhs.itsRealFundNotional;
	itsCallFreq					= (rhs.itsCallFreq);
	itsCallNotice				= (rhs.itsCallNotice);
	itsCallCal					= (rhs.itsCallCal);
	itsCallFees					= (rhs.itsCallFees);
	itsFundFreq					= (rhs.itsFundFreq);
	itsFundSpread				= (rhs.itsFundSpread);
	itsFundDayCount				= (rhs.itsFundDayCount);
	itsFundAdjRule				= (rhs.itsFundAdjRule);
	itsFundIntRule				= (rhs.itsFundIntRule);
	itsFundStubRule				= (rhs.itsFundStubRule);
	itsCpnSpread				= (rhs.itsCpnSpread);
	itsCpnPayFreq				= (rhs.itsCpnPayFreq);
	itsCpnResetCal				= (rhs.itsCpnResetCal);
	itsCpnPayCal				= (rhs.itsCpnPayCal);
	itsCpnAdjRule				= (rhs.itsCpnAdjRule);
	itsCpnIntRule				= (rhs.itsCpnIntRule);
	itsCpnStubRule				= (rhs.itsCpnStubRule);
	itsBoostedIndexType			= (rhs.itsBoostedIndexType);
	itsBoostedFixRate			= (rhs.itsBoostedFixRate);
	itsBoostedVarTerm			= (rhs.itsBoostedVarTerm);
	itsBoostedResetGap			= (rhs.itsBoostedResetGap);
	itsBoostedResetTiming		= (rhs.itsBoostedResetTiming);
	itsCpnResetGap				= (rhs.itsCpnResetGap);
	itsBoostedDayCount			= (rhs.itsBoostedDayCount);
	itsBoostedAdjRule			= (rhs.itsBoostedAdjRule);
	itsBoostedIntRule			= (rhs.itsBoostedIntRule);
	itsCpnBarrierDown			= (rhs.itsCpnBarrierDown);
	itsCpnBarrierUp				= (rhs.itsCpnBarrierUp);
	itsCpnResetFreq				= (rhs.itsCpnResetFreq);
	itsCpnResetTiming			= (rhs.itsCpnResetTiming);
	itsRefIndexType1			= (rhs.itsRefIndexType1);
	itsRefTerm1					= (rhs.itsRefTerm1);
	itsRefDayCount1				= (rhs.itsRefDayCount1);
	itsRefCoeff1				= (rhs.itsRefCoeff1);
	itsIsPortfolioNoticeDays	= (rhs.itsIsPortfolioNoticeDays);
	itsProductsToPrice			= (rhs.itsProductsToPrice);
	itsHasBeenPriced			= (rhs.itsHasBeenPriced);
	itsBermuda1Price			= (rhs.itsBermuda1Price);
	itsBermuda2Price			= (rhs.itsBermuda2Price);
	itsFundingPrice				= (rhs.itsFundingPrice);
	itsCorridorPrice			= (rhs.itsCorridorPrice);
	itsExoticSwapPrice			= (rhs.itsExoticSwapPrice);
	itsFwdFundingPrice			= (rhs.itsFwdFundingPrice);
	itsFwdCorridorPrice			= (rhs.itsFwdCorridorPrice);
	itsCorridorLegPrices		= (rhs.itsCorridorLegPrices);
	itsSwaption					= (rhs.itsSwaption ? (ARM_Swaption*) rhs.itsSwaption->Clone() : NULL);
	itsOptionPortfolio			= (rhs.itsOptionPortfolio ? (ARM_OptionPortfolio *) rhs.itsOptionPortfolio->Clone() : NULL);
	itsSigmaPf					= (rhs.itsSigmaPf);
	itsMrsPf					= (rhs.itsMrsPf);
	itsBetaPf					= (rhs.itsBetaPf);
	itsToCalibrateMeanReversion = (rhs.itsToCalibrateMeanReversion);
	itsToCalibrateBeta          = (rhs.itsToCalibrateBeta);
	itsReCalibMrs				= (rhs.itsReCalibMrs);
	itsReCalibBeta				= (rhs.itsReCalibBeta);
	itsFRMModel					= (rhs.itsFRMModel ? (ARM_FRMModel *) rhs.itsFRMModel->Clone() : NULL);
	itsVolType                  = (rhs.itsVolType);
	itsFlagToGenerateOSWATM     = (rhs.itsFlagToGenerateOSWATM);
	itsMeanRevMin               = (rhs.itsMeanRevMin);
	itsMeanRevMax               = (rhs.itsMeanRevMax);
	itsBetaMin                  = (rhs.itsBetaMin);
	itsBetaMax                  = (rhs.itsBetaMax);
	itsNbSteps                  = (rhs.itsNbSteps);
	itsMeanRev                  = (rhs.itsMeanRev);
	itsSigmaCurve               = (rhs.itsSigmaCurve ? (ARM_Vector *) rhs.itsSigmaCurve->Clone() : NULL);
	itsBetaOrShift              = (rhs.itsBetaOrShift ? (ARM_Vector *) rhs.itsBetaOrShift->Clone() : NULL);
	itsCalibSecPFParams			= (rhs.itsCalibSecPFParams ? (ARM_Vector *) rhs.itsCalibSecPFParams->Clone() : NULL);
	itsExerDateStrip			= (rhs.itsExerDateStrip);
	itsSFRMStdCalib				= (rhs.itsSFRMStdCalib);
}

////////////////////////////////////////////////////
///	Class  : ARM_CRACalculator
///	Routine: CleanUp
///	Returns: void
///	Action : delete pointors
////////////////////////////////////////////////////
void ARM_CRACalculator::CleanUp()
{
	if (itsSwaption)
	{
		delete itsSwaption;
		itsSwaption = NULL;
	}

	if (itsOptionPortfolio)
	{
		delete itsOptionPortfolio;
		itsOptionPortfolio = NULL;
	}

	if (itsCalibSecPFParams)
	{
		delete itsCalibSecPFParams;
		itsCalibSecPFParams = NULL;
	}

	if (itsSigmaCurve)
	{
		delete itsSigmaCurve;
		itsSigmaCurve = NULL;
	}

	if (itsBetaOrShift)
	{
		delete itsBetaOrShift;
		itsBetaOrShift = NULL;
	}

	if (itsFRMModel)
	{
		delete itsFRMModel;
		itsFRMModel = NULL;
	}
}

ARM_CRACalculator::ARM_CRACalculator(const ARM_MarketData_ManagerRep& mktDataManager)
                  :	ARM_GenCalculator(mktDataManager)
{
	itsSwaption = NULL;
	itsOptionPortfolio = NULL;
    itsCalibSecPFParams = NULL;

    itsFRMModel = NULL;

    itsSigmaCurve = NULL;
    itsBetaOrShift = NULL;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRACalculator 
///	Routine: Constructor 
///	Returns: void 
///	Action : builds the object 
/////////////////////////////////////////////////////////////////
ARM_CRACalculator::ARM_CRACalculator( ARM_Currency& ccy,
									  ARM_Date&	startDate,
									  ARM_Date&	endDate,
									  int payRec,
									  ARM_ReferenceValue&  notional,
					  				  int callFreq,
									  int callNotice,
									  string callCal,
									  ARM_ReferenceValue&  callFees,
									  int fundFreq,
									  ARM_ReferenceValue& fundSpread,
									  int fundDayCount,
									  ARM_ReferenceValue& cpnSpread,
									  int cpnPayFreq,
									  string cpnResetCal,
									  string cpnPayCal,
									  int boostedIndexType,
									  ARM_ReferenceValue& boostedFixRate,
									  string boostedVarTerm,
									  int boostedResetGap,
									  int boostedResetTiming,
									  int boostedDayCount,
									  int boostedAdjRule,
									  int boostedIntRule,
								 	  ARM_ReferenceValue& cpnBarrierDown,
									  ARM_ReferenceValue& cpnBarrierUp,
									  int cpnResetFreq,
									  int cpnResetTiming,
									  int cpnResetGap,
									  int refIndexType,
									  string refTerm,
									  int refDayCount,
									  double refCoeff,
									  double meanRevMin,
                                      double meanRevMax,
                                      double betaMin,
                                      double betaMax,
                                      ARM_Vector* calibSecPFParams,
                                      int nbSteps,
                                      int flagToGenerateOSWATM,
									  ARM_StringVector& mdmKeys,
									  const ARM_MarketData_ManagerRep& mktDataManager,
									  const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice,
									  int reCalibMrs,
									  int reCalibBeta,
									  bool isStdCalib,
									  ARM_Currency* fundccy,
									  const ARM_ReferenceValue& fundNominal )
: ARM_GenCalculator(mktDataManager),
itsCcy(ccy),
itsStartDate(startDate),
itsEndDate(endDate),
itsPayRec(payRec),
itsNotional(notional),
itsRealFundNotional(fundNominal),
itsCallFreq(callFreq),
itsCallNotice(callNotice),
itsCallCal(callCal),
itsCallFees(callFees),
itsFundFreq(fundFreq),
itsFundSpread(fundSpread),
itsFundDayCount(fundDayCount),
itsFundAdjRule(K_MOD_FOLLOWING),
itsFundIntRule(K_ADJUSTED),
itsFundStubRule(K_SHORTSTART),
itsCpnSpread(cpnSpread),
itsCpnPayFreq(cpnPayFreq),
itsCpnResetCal(cpnResetCal),
itsCpnPayCal(cpnPayCal),
itsCpnAdjRule(K_MOD_FOLLOWING),
itsCpnIntRule(K_ADJUSTED),
itsCpnStubRule(K_SHORTSTART),
itsBoostedIndexType(boostedIndexType),
itsBoostedFixRate(boostedFixRate),
itsBoostedVarTerm(boostedVarTerm),
itsBoostedResetGap(boostedResetGap),
itsBoostedResetTiming(boostedResetTiming),
itsBoostedDayCount(boostedDayCount), 
itsBoostedAdjRule(boostedAdjRule),
itsBoostedIntRule(boostedIntRule),
itsCpnBarrierDown(cpnBarrierDown),
itsCpnBarrierUp(cpnBarrierUp),
itsCpnResetFreq(cpnResetFreq),
itsCpnResetTiming(cpnResetTiming),
itsCpnResetGap(cpnResetGap),
itsRefIndexType1(refIndexType),
itsRefTerm1(refTerm),
itsRefDayCount1(refDayCount),
itsRefCoeff1(refCoeff),
itsIsPortfolioNoticeDays(false),
itsHasBeenPriced(false),
itsProductsToPrice(productsToPrice),
itsBermuda1Price(0.0),
itsBermuda2Price(0.0),
itsFundingPrice(0.0),
itsCorridorPrice(0.0),
itsFwdFundingPrice(0.0),
itsFwdCorridorPrice(0.0),
itsCorridorLegPrices(std::vector<double>(1,0.0)),
itsExoticSwapPrice(0.0),
itsSwaption(NULL),
itsOptionPortfolio(NULL),
itsSigmaPf(ARM_StdPortfolioPtr(NULL)),
itsMrsPf(ARM_StdPortfolioPtr(NULL)),
itsBetaPf(ARM_StdPortfolioPtr(NULL)),
itsToCalibrateMeanReversion(0),
itsToCalibrateBeta(0),
itsFRMModel(NULL),
itsVolType(K_DIAG),
itsMeanRevMin(meanRevMin),
itsMeanRevMax(meanRevMax),
itsBetaMin(betaMin),
itsBetaMax(betaMax),
itsReCalibMrs(reCalibMrs),
itsReCalibBeta(reCalibBeta),
itsCalibSecPFParams(calibSecPFParams ? (ARM_Vector *) calibSecPFParams->Clone() : NULL),
itsNbSteps(nbSteps),
itsMeanRev(0.0),
itsFlagToGenerateOSWATM(flagToGenerateOSWATM),
itsSigmaCurve(NULL),
itsBetaOrShift(NULL),
itsSFRMStdCalib(isStdCalib)
{
	if ( boostedIndexType == Fixed )
	{
		int period = (int)12/(int)ccy.GetFixedPayFreq();
		if ( period == 12 )
		{
			itsBoostedVarTerm = string("1Y");
		}
		else
		{
			char str[2];
			sprintf(str, "%d", period);
			itsBoostedVarTerm = string(str) + string("M");
		}
	}
	
	SetCurrencyUnit(&ccy); 
	SetFundingCcy(fundccy ? *fundccy: const_cast< ARM_Currency& >(ccy) );

	SetKeys(mdmKeys);

	CheckCRAInputs();

	CheckData();

	ARM_StringVector pricedColumns = PricedColumnNames();

	ARM_CstManagerPtr cstManagerPtr = CreateCstManager();

	Initialize();

    CreateAndSetDealDescriptionAndTimeIt("", pricedColumns, cstManagerPtr);

	CreateAndSetModel();

	CreateAndSetCalibration(); 
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRACalculator
///	Routine: Constructor 
///	Returns: void 
///	Action : Simple Contructor just to transfer 
/// datas useful for dealdescription
/////////////////////////////////////////////////////////////////
ARM_CRACalculator::ARM_CRACalculator(   const ARM_Currency& ccy,
										const ARM_Date&	startDate,
										const ARM_Date&	endDate,
										int payRec,
										const ARM_ReferenceValue&  notional,
										int callFreq,
					   				    int callNotice,
				   					    const string& callCal,
						  			    const ARM_ReferenceValue&  callFees,
						 		        int fundFreq,
						   			    const ARM_ReferenceValue& fundSpread,
										int fundDayCount,
					   				    const ARM_ReferenceValue& cpnSpread,
										int cpnPayFreq,
										const string& cpnResetCal,
										const string& cpnPayCal,
										int boostedIndexType,
										const ARM_ReferenceValue& boostedFixRate,
										const string& boostedVarTerm,
										int boostedResetGap,
										int boostedResetTiming,
										int boostedDayCount,
										int boostedAdjRule,
										int boostedIntRule,
										const ARM_ReferenceValue& cpnBarrierDown,
										const ARM_ReferenceValue& cpnBarrierUp,
										int cpnResetFreq,
										int cpnResetTiming,
										int cpnResetGap,
										int refIndexType1,
										const string& refTerm1,
									    int refDayCount1,
										double refCoeff1,
										bool  isPortfolioNoticeDays,
										const ARM_StringVector& mdmKeys,
										const ARM_MarketData_ManagerRep& mktDataManager,
										const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice,
										int reCalibMrs,
										int reCalibBeta,
										bool isStdCalib,
										ARM_Currency* fundccy,
										const ARM_ReferenceValue& fundNominal )
: ARM_GenCalculator(mktDataManager),
itsCcy(ccy),
itsStartDate(startDate),
itsEndDate(endDate),
itsPayRec(payRec),
itsNotional(notional),
itsRealFundNotional(fundNominal),
itsCallFreq(callFreq),
itsCallNotice(callNotice),
itsCallCal(callCal),
itsCallFees(callFees),
itsFundFreq(fundFreq),
itsFundSpread(fundSpread),
itsFundDayCount(fundDayCount),
itsFundAdjRule(K_MOD_FOLLOWING),
itsFundIntRule(K_ADJUSTED),
itsFundStubRule(K_SHORTSTART),
itsCpnSpread(cpnSpread),
itsCpnPayFreq(cpnPayFreq),
itsCpnResetCal(cpnResetCal),
itsCpnPayCal(cpnPayCal),
itsCpnAdjRule(K_MOD_FOLLOWING),
itsCpnIntRule(K_ADJUSTED),
itsCpnStubRule(K_SHORTSTART),
itsBoostedIndexType(boostedIndexType),
itsBoostedFixRate(boostedFixRate),
itsBoostedVarTerm(boostedVarTerm),
itsBoostedResetGap(boostedResetGap),
itsBoostedResetTiming(boostedResetTiming),
itsBoostedDayCount(boostedDayCount), 
itsBoostedAdjRule(boostedAdjRule),
itsBoostedIntRule(boostedIntRule),
itsCpnBarrierDown(cpnBarrierDown),
itsCpnBarrierUp(cpnBarrierUp),
itsCpnResetFreq(cpnResetFreq),
itsCpnResetTiming(cpnResetTiming),
itsCpnResetGap(cpnResetGap),
itsRefIndexType1(refIndexType1),
itsRefTerm1(refTerm1),
itsRefDayCount1(refDayCount1),
itsRefCoeff1(refCoeff1),
itsIsPortfolioNoticeDays(isPortfolioNoticeDays),
itsHasBeenPriced(false),
itsProductsToPrice(productsToPrice),
itsBermuda1Price(0.0),
itsBermuda2Price(0.0),
itsFundingPrice(0.0),
itsCorridorPrice(0.0),
itsFwdFundingPrice(0.0),
itsFwdCorridorPrice(0.0),
itsCorridorLegPrices(std::vector<double>(1,0.0)),
itsExoticSwapPrice(0.0),
itsSwaption(NULL),
itsOptionPortfolio(NULL),
itsSigmaPf(ARM_StdPortfolioPtr(NULL)),
itsMrsPf(ARM_StdPortfolioPtr(NULL)),
itsBetaPf(ARM_StdPortfolioPtr(NULL)),
itsToCalibrateMeanReversion(0),
itsToCalibrateBeta(0),
itsFRMModel(NULL),
itsVolType(K_DIAG),
itsMeanRevMin(0.0),
itsMeanRevMax(0.0),
itsBetaMin(0.0),
itsBetaMax(0.0),
itsReCalibMrs(reCalibMrs),
itsReCalibBeta(reCalibBeta),
itsCalibSecPFParams(NULL),
itsNbSteps(1),
itsMeanRev(0.0),
itsFlagToGenerateOSWATM(true),
itsSigmaCurve(NULL),
itsBetaOrShift(NULL),
itsSFRMStdCalib(isStdCalib)
{
	SetCurrencyUnit(const_cast<ARM_Currency*>(&ccy)); 
	SetFundingCcy(fundccy ? *fundccy: const_cast< ARM_Currency& >(ccy) );

	SetKeys(mdmKeys);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRACalculator 
///	Routine: Summit Constructor 
///	Returns: void 
///	Action : builds the object 
/////////////////////////////////////////////////////////////////
ARM_CRACalculator::ARM_CRACalculator( ARM_OptionPortfolio* optPortfolio,
						              double meanRevMin,
                                      double meanRevMax,
                                      double betaMin,
                                      double betaMax,
                                      ARM_Vector* calibSecPFParams,
                                      int nbSteps,
                                      int flagToGenerateOSWATM,
									  ARM_StringVector& mdmKeys,
									  const ARM_MarketData_ManagerRep& mktDataManager,
									  const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice,
									  int reCalibMrs,
									  int reCalibBeta,
									  bool isStdCalib)   							  
: ARM_GenCalculator(mktDataManager),
itsHasBeenPriced(false),
itsProductsToPrice(productsToPrice),
itsBermuda1Price(0.0),
itsBermuda2Price(0.0),
itsFundingPrice(0.0),
itsCorridorPrice(0.0),
itsFwdFundingPrice(0.0),
itsFwdCorridorPrice(0.0),
itsCorridorLegPrices(std::vector<double>(1,0.0)),
itsExoticSwapPrice(0.0),
itsSwaption(NULL),
itsOptionPortfolio(NULL),
itsSigmaPf(ARM_StdPortfolioPtr(NULL)),
itsMrsPf(ARM_StdPortfolioPtr(NULL)),
itsBetaPf(ARM_StdPortfolioPtr(NULL)),
itsToCalibrateMeanReversion(0),
itsToCalibrateBeta(0),
itsFRMModel(NULL),
itsVolType(K_DIAG),
itsMeanRevMin(meanRevMin),
itsMeanRevMax(meanRevMax),
itsBetaMin(betaMin),
itsBetaMax(betaMax),
itsReCalibMrs(reCalibMrs),
itsReCalibBeta(reCalibBeta),
itsCalibSecPFParams(calibSecPFParams ? (ARM_Vector *) calibSecPFParams->Clone() : NULL),
itsNbSteps(nbSteps),
itsMeanRev(0.0),
itsFlagToGenerateOSWATM(flagToGenerateOSWATM),
itsSigmaCurve(NULL),
itsBetaOrShift(NULL),
itsIsPortfolioNoticeDays(true),
itsSFRMStdCalib(isStdCalib)
{
	itsCcy = *(optPortfolio->GetLiborLeg()->GetCurrencyUnit());
	
	SetCurrencyUnit(&itsCcy); 

	SetKeys(mdmKeys);

	CheckCRAInputs();

	CheckData();

	GenerateProductDescription((ARM_OptionPortfolio *) optPortfolio);

	Initialize();

	ARM_StringVector pricedColumns = PricedColumnNames();

	ARM_CstManagerPtr cstManagerPtr = CreateCstManager();

    CreateAndSetDealDescriptionAndTimeIt("", pricedColumns, cstManagerPtr);

	CreateAndSetModel();

    CreateAndSetCalibration(); 
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRACalculator 
///	Routine: Summit Constructor 
///	Returns: void 
///	Action : builds a calculator without market data 
/////////////////////////////////////////////////////////////////
ARM_CRACalculator::ARM_CRACalculator( const ARM_Date& asOfDate,
									  ARM_Security* security)
: ARM_GenCalculator(asOfDate),
itsCcy(*security->GetCurrencyUnit()),
itsHasBeenPriced(false),
itsProductsToPrice(NULL),
itsBermuda1Price(0.0),
itsBermuda2Price(0.0),
itsFundingPrice(0.0),
itsCorridorPrice(0.0),
itsFwdFundingPrice(0.0),
itsFwdCorridorPrice(0.0),
itsCorridorLegPrices(std::vector<double>(1,0.0)),
itsExoticSwapPrice(0.0),
itsSwaption(NULL),
itsOptionPortfolio(NULL),
itsSigmaPf(ARM_StdPortfolioPtr(NULL)),
itsMrsPf(ARM_StdPortfolioPtr(NULL)),
itsBetaPf(ARM_StdPortfolioPtr(NULL)),
itsToCalibrateMeanReversion(0),
itsToCalibrateBeta(0),
itsFRMModel(NULL),
itsVolType(K_DIAG),
itsMeanRevMin(0.0),
itsMeanRevMax(0.0),
itsBetaMin(0.0),
itsBetaMax(0.0),
itsReCalibMrs(-1),
itsReCalibBeta(-1),
itsCalibSecPFParams(NULL),
itsNbSteps(0),
itsMeanRev(0.0),
itsFlagToGenerateOSWATM(1),
itsSigmaCurve(NULL),
itsBetaOrShift(NULL),
itsIsPortfolioNoticeDays(true),
itsSFRMStdCalib(true)
{
	SetCurrencyUnit(security->GetCurrencyUnit()); 
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRACalculator
///	Routine: Copy constructor
///	Returns: a built object
///	Action : builds the object
/////////////////////////////////////////////////////////////////
ARM_CRACalculator::ARM_CRACalculator(const ARM_CRACalculator& rhs)
:	ARM_GenCalculator(rhs)
{
	CopyNoCleanUp(rhs);
}

void ARM_CRACalculator::Set( ARM_Currency& ccy,
							  ARM_Date&	startDate,
							  ARM_Date&	endDate,
							  int payRec,
							  ARM_ReferenceValue&  notional,
					  		  int callFreq,
							  int callNotice,
							  string callCal,
							  ARM_ReferenceValue&  callFees,
							  int fundFreq,
							  ARM_ReferenceValue& fundSpread,
							  int fundDayCount,
							  ARM_ReferenceValue& cpnSpread,
							  int cpnPayFreq,
							  string cpnResetCal,
							  string cpnPayCal,
							  int boostedIndexType,
							  ARM_ReferenceValue& boostedFixRate,
							  string boostedVarTerm,
							  int boostedResetGap,
							  int boostedResetTiming,
							  int boostedDayCount,
							  int boostedAdjRule,
							  int boostedIntRule,
							  ARM_ReferenceValue& cpnBarrierDown,
							  ARM_ReferenceValue& cpnBarrierUp,
							  int cpnResetFreq,
							  int cpnResetTiming,
							  int cpnResetGap,
							  int refIndexType1,
							  string refTerm1,
							  int refDayCount1,
							  double refCoeff1,
                              double meanRevMin,
                              double meanRevMax,
                              double betaMin,
                              double betaMax,
                              int nbSteps,
                              ARM_Vector* calibSecPFParams,
							  int flagToGenerateOSWATM,
							  ARM_StringVector& mdmKeys,
							  const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice,
							  int reCalibMrs,
							  int reCalibBeta,
							  int volType,
							  bool isStdCalib)
{
	SetCurrencyUnit(&ccy); 
	SetCcy(ccy);
	SetStartDate(startDate);
	SetEndDate(endDate);
	SetPayRec(payRec);
	SetNotional(notional);
	SetCallFreq(callFreq);
	SetCallNotice(callNotice);
	SetCallCal(callCal);
	SetCallFees(callFees);
	SetFundFreq(fundFreq);
	SetFundSpread(fundSpread);
	SetFundDayCount(fundDayCount);
	SetCpnSpread(cpnSpread);
	SetCpnPayFreq(cpnPayFreq);
	SetCpnResetCal(cpnResetCal);
	SetCpnPayCal(cpnPayCal);
	SetBoostedIndexType(boostedIndexType);
	SetBoostedFixRate(boostedFixRate);
	SetBoostedVarTerm(boostedVarTerm);
	SetBoostedResetGap(boostedResetGap);
	SetBoostedResetTiming(boostedResetTiming);
	SetBoostedDayCount(boostedDayCount);
	SetBoostedAdjRule(boostedAdjRule);
	SetBoostedIntRule(boostedIntRule);
	SetCpnBarrierDown(cpnBarrierDown);
	SetCpnBarrierUp(cpnBarrierUp);
	SetCpnResetFreq(cpnResetFreq);
	SetCpnResetTiming(cpnResetTiming);
	SetCpnResetGap(cpnResetGap);
	SetRefIndexType1(refIndexType1);
	SetRefTerm1(refTerm1);
	SetRefDayCount1(refDayCount1);
	SetRefCoeff1(refCoeff1);
	SetMeanRevMin(meanRevMin);
	SetMeanRevMax(meanRevMax);
	SetBetaMin(betaMin);
	SetBetaMax(betaMax);
	SetCalibSecPFParams(calibSecPFParams ? (ARM_Vector*) calibSecPFParams->Clone() : NULL);
	SetNbSteps(nbSteps);
	SetFlagToGenerateOSWATM(flagToGenerateOSWATM);
	SetKeys(mdmKeys);
	SetProductsToPrice(productsToPrice);
	SetReCalibMrs(reCalibMrs);
	SetReCalibBeta(reCalibBeta);
	SetVolType(volType);
	SetSfrmStdCalib(isStdCalib);
	itsMeanRev = 0.0;
	itsIsPortfolioNoticeDays = false;
	itsHasBeenPriced = false;
	itsBermuda1Price = 0.0;
	itsBermuda2Price = 0.0;
	itsFundingPrice = 0.0;
	itsCorridorPrice = 0.0;
	itsFwdFundingPrice = 0.0;
	itsFwdCorridorPrice = 0.0;
	itsCorridorLegPrices = std::vector<double>(1,0.0);
	itsExoticSwapPrice = 0.0;
	itsSwaption = NULL;
	itsOptionPortfolio = NULL;
	itsSigmaPf = ARM_StdPortfolioPtr(NULL);
	itsMrsPf = ARM_StdPortfolioPtr(NULL);
	itsBetaPf = ARM_StdPortfolioPtr(NULL);
	itsToCalibrateMeanReversion = 0;
	itsToCalibrateBeta = 0;
	itsReCalibMrs = reCalibMrs;
	itsReCalibBeta = reCalibBeta;
	itsFRMModel = NULL;
	itsSigmaCurve = NULL;
	itsBetaOrShift = NULL;
}

void ARM_CRACalculator::Set(ARM_OptionPortfolio* optPortfolio,
							double meanRevMin,
							double meanRevMax,
							double betaMin,
							double betaMax,
							ARM_Vector* calibSecPFParams,
							int nbSteps,
							int flagToGenerateOSWATM,
							ARM_StringVector& mdmKeys,
							const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice,
							int reCalibMrs,
							int reCalibBeta,
							int volType,
							bool isStdCalib)
{
	SetCurrencyUnit(optPortfolio->GetLiborLeg()->GetCurrencyUnit()); 
	SetCcy(*(optPortfolio->GetLiborLeg()->GetCurrencyUnit()));
	SetOptionPortfolio(optPortfolio);
	SetMeanRevMin(meanRevMin);
	SetMeanRevMax(meanRevMax);
	SetBetaMin(betaMin);
	SetBetaMax(betaMax);
	SetCalibSecPFParams(calibSecPFParams ? (ARM_Vector*) calibSecPFParams->Clone() : NULL);
	SetNbSteps(nbSteps);
	SetFlagToGenerateOSWATM(flagToGenerateOSWATM);
	SetKeys(mdmKeys);
	SetProductsToPrice(productsToPrice);
	SetReCalibMrs(reCalibMrs);
	SetReCalibBeta(reCalibBeta);
	SetVolType(volType);
	SetSfrmStdCalib(isStdCalib);
	itsMeanRev = 0.0;
	itsHasBeenPriced = false;
	itsBermuda1Price = 0.0;
	itsBermuda2Price = 0.0;
	itsFundingPrice = 0.0;
	itsCorridorPrice = 0.0;
	itsFwdFundingPrice = 0.0;
	itsFwdCorridorPrice = 0.0;
	itsCorridorLegPrices = std::vector<double>(1,0.0);
	itsExoticSwapPrice = 0.0;
	itsSigmaPf = ARM_StdPortfolioPtr(NULL);
	itsMrsPf = ARM_StdPortfolioPtr(NULL);
	itsBetaPf = ARM_StdPortfolioPtr(NULL);
	itsToCalibrateMeanReversion = 0;
	itsToCalibrateBeta = 0;
	itsFRMModel = NULL;
	itsReCalibMrs = reCalibMrs;
	itsReCalibBeta = reCalibBeta;
	itsSigmaCurve = NULL;
	itsBetaOrShift = NULL;
	itsIsPortfolioNoticeDays = true;
	itsSwaption = NULL;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRACalculator
///	Routine: Destructor
///	Returns: void
///	Action : destroys the object
/////////////////////////////////////////////////////////////////
ARM_CRACalculator::~ARM_CRACalculator()
{
	CleanUp();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRACalculator
///	Routine: assignment operator
///	Returns: 
///	Action : 
/////////////////////////////////////////////////////////////////
ARM_CRACalculator& ARM_CRACalculator::operator=(const ARM_CRACalculator& rhs)
{
	if( this != & rhs )
	{
		ARM_GenCalculator::operator=(rhs);
        CleanUp();
        CopyNoCleanUp(rhs);
	}
	return *this;
}

////////////////////////////////////////////////////
///	Class   : ARM_CRACalculator
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
///           Call copy constructor
////////////////////////////////////////////////////
ARM_Object* ARM_CRACalculator::Clone() const
{
	return(new ARM_CRACalculator(*this));
}

////////////////////////////////////////////////////
///    Class  : ARM_CRACalculator
///    Routine: Pre_InitialiseFRMModel
///    Returns: 
///    Action : To Preinitialize the ARM_FRMModel
////////////////////////////////////////////////////
void ARM_CRACalculator::Pre_InitialiseFRMModel()
{
	ARM_Model* capModel = dynamic_cast< ARM_Model* >( GetMktDataManager()->GetData(GetKeys()[CfModelKey]) );
    if (dynamic_cast<ARM_BSModel*>(capModel))
	  capModel = dynamic_cast<ARM_BSModel*>(capModel);
	else if (dynamic_cast<ARM_BSSmiledModel*>(capModel))
	  capModel = dynamic_cast<ARM_BSSmiledModel*>(capModel);

	vector<double> calibParams(5);
	calibParams[0] = 0.0;
	calibParams[1] = itsMeanRevMin;
	calibParams[2] = itsMeanRevMax;
	calibParams[3] = itsBetaMin;
	calibParams[4] = itsBetaMax;

	int size = calibParams.size();
    ARM_Vector temparam(size-1);
    for (int i = 1; i < size; i++)
        temparam[i-1] = (calibParams)[i];

    long nbFactors = 1; // 1 Factor So far!

    ARM_CapFloor* capFloor = (ARM_CapFloor *) itsSigmaPf->GetAsset(0);
    ARM_INDEX_TYPE liborType = capFloor->GetSwapLeg()->GetIRIndex()->GetIndexType();

    ARM_Vector* clonedSigmaCrv = itsSigmaCurve ? (ARM_Vector *) itsSigmaCurve->Clone() : NULL;

	int szBeta = itsBetaOrShift->GetSize();
	double cstbeta = (itsBetaMax + itsBetaMin)/2.0; 
	ARM_Vector* clonedBetaOrShift = NULL;

	if (!itsReCalibBeta)
		clonedBetaOrShift = itsBetaOrShift ? (ARM_Vector *) itsBetaOrShift->Clone(): NULL;
	else
		clonedBetaOrShift = new ARM_Vector(szBeta, cstbeta);
    
    ARM_Date asOfDate = capModel->GetStartDate();
    ARM_FRMHWVol* FRMHWVol = new ARM_FRMHWVol(asOfDate,
                                              liborType,
                                              (ARM_Portfolio*)(&*itsSigmaPf),
                                              (ARM_Portfolio*)(&*itsBetaPf),
                                              clonedSigmaCrv,
                                              clonedBetaOrShift,
                                              0,
                                              (itsMeanRevMin + itsMeanRevMax) / 2.0,
                                              nbFactors,
                                              NULL,
                                              itsVolType);

    if (itsReCalibBeta)
       itsFRMModel = new ARM_FRMModel(asOfDate, 
                                      capModel->GetZeroCurve(),
                                      FRMHWVol,
                                      (ARM_Portfolio*)(&*itsSigmaPf),
                                      (ARM_Portfolio*)(&*itsMrsPf),
                                      (ARM_Portfolio*)(&*itsBetaPf),
                                      &temparam,
                                      0,
                                      capModel,
                                      1);
    else
       itsFRMModel = new ARM_FRMModel(asOfDate, 
                                      capModel->GetZeroCurve(),
                                      FRMHWVol,
                                      (ARM_Portfolio*)(&*itsSigmaPf),
                                      (ARM_Portfolio*)(&*itsMrsPf),
                                      (ARM_Portfolio*)(&*itsBetaPf),
                                      &temparam,
                                      0);
}

////////////////////////////////////////////////////
///	Class   : ARM_CRACalculator
///	Routines: Initialize
///	Returns : void
///	Action  : Create calibration portfolios
///           Compute Sigma and Beta
////////////////////////////////////////////////////
void ARM_CRACalculator::Initialize()
{
	try
	{	
		if (!itsOptionPortfolio)
		{
			itsOptionPortfolio = GenerateOptionPortfolio();
		}

		itsOptionPortfolio->SetFlagToGenerateOSWATM(itsFlagToGenerateOSWATM ? "YES" : "NO");
		itsOptionPortfolio->SetPFCalibParams(itsCalibSecPFParams);
	
		ARM_Model* capModel = dynamic_cast< ARM_Model* >( GetMktDataManager()->GetData(GetKeys()[CfModelKey]) );
		if (dynamic_cast<ARM_BSModel*>(capModel))
		  capModel = dynamic_cast<ARM_BSModel*>(capModel);
		else if (dynamic_cast<ARM_BSSmiledModel*>(capModel))
		  capModel = dynamic_cast<ARM_BSSmiledModel*>(capModel);

		ARM_Date asofdate = capModel->GetStartDate();

		double meanRevMin = itsMeanRevMin;
		double meanRevMax = itsMeanRevMax;

		itsToCalibrateMeanReversion = fabs(meanRevMax - meanRevMin) > K_NEW_DOUBLE_TOL ? true : false;
		itsReCalibMrs = itsReCalibMrs == -1 ? itsToCalibrateMeanReversion : itsReCalibMrs;

		if(itsReCalibMrs && (!itsToCalibrateMeanReversion))
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
							 "calibration flag are not valid : minMRS and MaxMrs should not be equal");

		itsToCalibrateMeanReversion = itsReCalibMrs;


		double betaOrShiftMin = itsBetaMin;
		double betaOrShiftMax = itsBetaMax;

		itsToCalibrateBeta = fabs(betaOrShiftMax - betaOrShiftMin) > K_NEW_DOUBLE_TOL ? true : false;
		itsReCalibBeta = itsReCalibBeta == -1 ? itsToCalibrateBeta : itsReCalibBeta;

		if(itsReCalibBeta && (!itsToCalibrateBeta))
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
								 "calibration flag are not valid : minBeta and MaxBeta should not be equal");

		itsToCalibrateBeta = itsReCalibBeta;

		long upDownBarrierFlag = itsOptionPortfolio->GetUpDownBarrierFlag();
		if((itsToCalibrateBeta) && (upDownBarrierFlag < 0))
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
								 "We don't calibrate corridors for ATM strikes ");

		//Sigma portfolio
		if (itsSigmaPf.IsNull())
		{
			itsSigmaPf = ARM_StdPortfolioPtr(itsOptionPortfolio->GenerateSigDiagCalibPF(asofdate,
											 itsOptionPortfolio->GetCapNbMonths(),
											 itsOptionPortfolio->GetPFStyle(),
											 itsOptionPortfolio->GetUpDownBarrierFlag()));

			size_t size = itsSigmaPf->size();

			for (int i = 0; i < size; ++i)
			{
				itsSigmaPf->GetAsset(i)->SetModelVariable(NULL);
				itsSigmaPf->GetAsset(i)->SetModel(capModel);

				double mktPrice = itsSigmaPf->GetAsset(i)->ComputePrice();
				double vega = itsSigmaPf->GetAsset(i)->ComputeSensitivity(K_VEGA);
				double curWeight = ( vega > 1.0e-5 ) ? 1.0 : 0.0;       

				itsSigmaPf->SetWeight(curWeight, i);
				itsSigmaPf->SetPrecision(vega, i);
				itsSigmaPf->SetPrice(mktPrice, i);
			}
		}

		//Beta portfolio
		if (itsToCalibrateBeta && itsBetaPf.IsNull())
		{
			itsBetaPf = ARM_StdPortfolioPtr(itsOptionPortfolio->GenerateBetaCalibPF(itsOptionPortfolio->GetSwitchFreqIfDaily(),
																					itsOptionPortfolio->GetUpDownBarrierFlag()));

			size_t size = itsBetaPf->size();
			for (int i = 0; i < size; ++i)
			{
				itsBetaPf->GetAsset(i)->SetModelVariable(NULL);
				itsBetaPf->GetAsset(i)->SetModel(capModel);

				double mktPrice = itsBetaPf->GetAsset(i)->ComputePrice();
				double vega = itsBetaPf->GetAsset(i)->ComputeSensitivity(K_VEGA);
				double curWeight = ( i == 0 ) ? 100 : 1.0;

				itsBetaPf->SetWeight(curWeight, i);
				itsBetaPf->SetPrecision(vega, i);
				itsBetaPf->SetPrice(mktPrice, i);
			}
		}
   
		//Mean reversion portfolio
		if (itsToCalibrateMeanReversion && itsMrsPf.IsNull())
		{
			int freqChngeIfDaily = 0;  // let it daily 

			itsMrsPf = ARM_StdPortfolioPtr(itsOptionPortfolio->GenerateMeanRevCalibPF(capModel,
													   itsOptionPortfolio->GetSwoptNbMonths(),
													   freqChngeIfDaily));

  			ARM_Model* swptModel = (ARM_Model*) (GetMktDataManager()->GetData(GetKeys()[OswModelKey]));
        
			size_t size = itsMrsPf->size();
			for (int i = 0; i < size; ++i)
			{
				itsMrsPf->GetAsset(i)->SetModelVariable(NULL);
				itsMrsPf->GetAsset(i)->SetModel(swptModel);

				double mktPrice = itsMrsPf->GetAsset(i)->ComputePrice();
				double vega = itsMrsPf->GetAsset(i)->ComputeSensitivity(K_VEGA);
				double curWeight = ( vega > 1.0e-5 ) ? 1.0 : 0.0;       

				itsMrsPf->SetWeight(curWeight, i);
				itsMrsPf->SetPrecision(vega, i);
				itsMrsPf->SetPrice(mktPrice, i);
			}
		}

		int sz = itsSigmaPf->GetSize();
		itsSigmaCurve = new ARM_Vector(sz); 

		int i;
		for (i = 0; i < sz; i++)
		{
			ARM_CapFloor* capFloor = (ARM_CapFloor *) itsSigmaPf->GetAsset(i);

			if (itsSigmaPf->GetWeights()->Elt(i))
			{
			   itsSigmaCurve->Elt(i) = capFloor->GetCapletAdjVols()->Elt(0)*100.0;
			}
			else
			{
			   double asofdate = capFloor->GetModel()->GetStartDate().GetJulian();
			   double optMat = (capFloor->GetSwapLeg()->GetResetDates()->Elt(0) - asofdate)/K_YEAR_LEN;
			   ARM_IRIndex* irIndex = capFloor->GetSwapLeg()->GetIRIndex();
			   int dayCount = irIndex->GetDayCount();

			   double optLetMat = CountYears(dayCount, capFloor->GetSwapLeg()->GetFwdRateStartDates()->Elt(0), 
											 capFloor->GetSwapLeg()->GetFwdRateEndDates()->Elt(0));
          
			   ARM_Model* model = capFloor->GetModel();

			   double vol;

			   if (dynamic_cast<ARM_BSModel*>(capModel))
				  vol = ((ARM_BSModel*)capModel)->GetVolatility()->ComputeVolatility(optMat, optLetMat);
			   else if (dynamic_cast<ARM_BSSmiledModel*>(capModel))
				  vol = ((ARM_BSSmiledModel*)capModel)->GetVolatility()->ComputeVolatility(optMat, optLetMat);
			   else
				  throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
								 "No model valid to compute volatility if ARM_Calibrator SFRM");

			   itsSigmaCurve->Elt(i) = vol;
			}
		}

		int szBeta = itsToCalibrateBeta ? itsBetaPf->GetSize() : itsSigmaPf->GetSize();
		double cstbeta = (itsBetaMax + itsBetaMin)/2.0; 
		itsBetaOrShift = new ARM_Vector(szBeta, cstbeta);
	}
	catch (Exception& x)
	{
	    x.DebugPrint();

		throw x;
	}
	catch (...)
	{
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in CraCalculator::Initialize" );
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRACalculator
///	Routine: CreateCstManager
///	Returns: void
///	Action : create the const manager (static data).
/////////////////////////////////////////////////////////////////
ARM_CstManagerPtr ARM_CRACalculator::CreateCstManager()
{
	try
	{
		double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

		vector<string>	cstNames; 
		cstNames.push_back("BoostedFixRateCst");
		cstNames.push_back("CpnSpreadCst");
		cstNames.push_back("BarrierDownCst");
		cstNames.push_back("BarrierUpCst");
		cstNames.push_back("NotionalCst");
		cstNames.push_back("FundSpreadCst");

		vector<ARM_GramFctorArg> cstVector;
		
		//BoostedFixRateCst
		ARM_Curve* tmpBoosted;
		tmpBoosted = RefValueToCurve(GetBoostedFixRate(), asOfDate);
		cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpBoosted))));

		if (itsBoostedIndexType != K_FIXED)
		{
			//CpnSpreadCst (Stepup left)
			ARM_Curve* tmpCpnSpread;
			tmpCpnSpread = RefValueToCurve(GetCpnSpread(), asOfDate);
			cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpCpnSpread))));
		}
		else
		{
			std::vector<double> dates(1, 0.0);
			std::vector<double> values(1, 0.0);
			ARM_Curve* tmpCpnSpread = new ARM_Curve(dates, values,new ARM_LinInterpCstExtrapolDble);
			cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpCpnSpread))));
		}

		//BarrierDownCst (Stepup left)
		ARM_Curve* tmpBarrierDown = RefValueToCurve(itsCpnBarrierDown, asOfDate);
		cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpBarrierDown))));

		//BarrierUpCst (Stepup left)
		ARM_Curve* tmpBarrierUp = RefValueToCurve(itsCpnBarrierUp, asOfDate);
		cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpBarrierUp))));

		//NotionalCst (Stepup right)
		ARM_Curve* tmpNotional = RefValueToCurve(itsNotional, asOfDate);
		cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpNotional))));

		//FundSpreadCst (Stepup left)
		ARM_Curve* tmpFundSpread = RefValueToCurve(itsFundSpread, asOfDate);
		cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpFundSpread))));

		ARM_CstManagerPtr cstManagerPtr = ARM_CstManagerPtr(new ARM_CstManager(cstNames, cstVector));

		return cstManagerPtr;
	}
	catch (Exception& x)
	{
	    x.DebugPrint();

		throw x;
	}
	catch (...)
	{
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in CraCalculator::CreateCstManager" );
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRACalculator
///	Routine: GetIndexType
///	Returns: the index type
///	Action : get the ARM index type of the coupon index
/////////////////////////////////////////////////////////////////
ARM_INDEX_TYPE ARM_CRACalculator::GetIndexType(int type)
{
	string varIndexTerm;

    if (type == Boosted)
		varIndexTerm = itsBoostedVarTerm;
	else if (type == RefIndex)
		varIndexTerm = itsRefTerm1;
		
    string liborTypeName(string(GetCcy().GetCcyName()) == string("EUR") ? "EURIBOR" : "LIBOR");
    if(varIndexTerm=="12M")
        varIndexTerm="1Y"; // because we cant convert EURIBOR1Y and not EURIBOR12M
    liborTypeName += varIndexTerm;
    return static_cast< ARM_INDEX_TYPE > (ARM_ArgConv_IndexType.GetNumber(liborTypeName));
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRACalculator
///	Routine: CheckCRAInputs
///	Returns: void
///	Action : check if Callable Range Accrual inputs are consistent
/////////////////////////////////////////////////////////////////
void ARM_CRACalculator::CheckCRAInputs()
{	
	//to fullfill when necessary.
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRACalculator
///	Routine: CheckMktData & CheckData
///	Returns: void
///	Action : check if data are consistent
/////////////////////////////////////////////////////////////////
void ARM_CRACalculator::CheckMktData()
{
    /// MdM datas checking
	ARM_ZeroCurve* ycCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));
    if(!ycCurve)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Yield Curve for key=" + GetKeys()[YcKey] + " is expected in the Market Data Manager");
    string Ccy(ycCurve->GetCurrencyUnit()->GetCcyName());

    ARM_BSModel* oswBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswModelKey]) );
    if(!oswBSModel)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : swaption B&S model for key=" + GetKeys()[OswModelKey] + " is expected in the Market Data Manager");

	ARM_BSModel* cfBSModel = dynamic_cast< ARM_BSModel* >(GetMktDataManager()->GetData(GetKeys()[CfModelKey]) );
    if(!cfBSModel)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : cap/floor B&S smiled model for key=" + GetKeys()[CfModelKey] + " is expected in the Market Data Manager");

	ARM_ModelParam* mrsParam = dynamic_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[MrsKey]));
    if(!mrsParam || mrsParam->GetType() != ARM_ModelParamType::MeanReversion)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : MRS Param for key=" + GetKeys()[MrsKey] + " is expected in the Market Data Manager");
}

void ARM_CRACalculator::CheckData()
{
	CheckMktData();
}

//////////////////////////////////////////////////////////////////////////////
///	Class  : ARM_CRACalculator
///	Routine: PricedColumnNames
///	Returns: ARM_StringVector
///	Action : create the priced column names of the deal description
//////////////////////////////////////////////////////////////////////////////
ARM_StringVector ARM_CRACalculator::PricedColumnNames() const
{
	ARM_StringVector pricedColumns;

	// effective nb of columns to price :
	int colSize = MIN(NbProductsToPrice, itsProductsToPrice.size());
	
	for (int i=0; i<colSize; i++)
	{
		if (itsProductsToPrice[i])
			pricedColumns.push_back(CRAColNamesTable[CRAProductToPriceColumns [i]]);
	}

	return pricedColumns;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRACalculator
///	Routine: ColumnNames
///	Returns: ARM_RowInfo
///	Action : create the column names of the deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_CRACalculator::ColumnNames() const
{
	// Number of Columns
	size_t colNamesSize;

	colNamesSize = sizeof(CRAColNamesTable)/sizeof(CRAColNamesTable[0]);

    vector< string > colNamesVec(colNamesSize);
    vector< ARM_GP_VALUE_TYPE > colTypeVec(colNamesSize, ARM_STRING); 

	for (size_t i = 0; i < colNamesSize; ++i)
	{
		colNamesVec[i] = CRAColNamesTable[i];
	}

    ARM_RowInfo rowInfo(colNamesVec,colTypeVec);

    return rowInfo;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRACalculator
///	Routine: MiddleRows
///	Returns: ARM_RowInfo
///	Action : create a row of a deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_CRACalculator::MiddleRows(size_t eventIdx, 
                                          const ARM_DateStripCombiner& datesStructure) const
{
	try
	{
		//Useful datas
		char* ccy	= GetCcy().GetCcyName();
		double maturityDate = GetEndDate().GetJulian();
		int intPayRec = GetPayRec();
		string payRec; 
		string invPayRec;
		if ( intPayRec == -1 )
		{
			payRec = "P";
			invPayRec = "R";
		}
		else
		{
			payRec = "R";
			invPayRec = "P";
		}
		int intFundFreq = GetFundFreq();
		string fundFreq = ARM_ArgConvReverse_MatFrequency.GetString(intFundFreq);
		int intFundDayCount = GetFundDayCount();
		string fundDayCount = ARM_ArgConvReverse_DayCount.GetString(intFundDayCount);
		int intBoosIndexType = GetBoostedIndexType();
		string boostedIndexType;
		if (intBoosIndexType == Fixed)
			boostedIndexType = "FIXED";
		else if (intBoosIndexType == Libor)
			boostedIndexType = "LIBOR";
		else if (intBoosIndexType == Cms)
			boostedIndexType = "CMS";
		int intCpnPayFreq = GetCpnPayFreq();
		string cpnPayFreq = ARM_ArgConvReverse_StdFrequency.GetString(intCpnPayFreq);
		string boostedVarTerm = GetBoostedVarIndexTerm();
		int intBoosResetTiming = GetBoostedResetTiming();
		string boostedResetTiming = ARM_ArgConvReverse_Timing.GetString(intBoosResetTiming);
		int intBoosPayDayCount = GetBoostedDayCount(); 
		string boostedPayDayCount = ARM_ArgConvReverse_DayCount.GetString(intBoosPayDayCount);
		int cpnResetGap = GetCpnResetGap();
		int intCpnIntRule = GetBoostedIntRule();
		string cpnIntRule = ARM_ArgConvReverse_IntRule.GetString(intCpnIntRule);
		int intCpnResetFreq = GetCpnResetFreq();
		string cpnResetFreq = ARM_ArgConvReverse_StdFrequency.GetString(intCpnResetFreq);
		int intCpnResetTiming = GetCpnResetTiming();
		string cpnResetTiming = ARM_ArgConvReverse_Timing.GetString(intCpnResetTiming);
		int intCpnIndexType1 = GetRefIndexType1();
		string cpnIndexType1 = ARM_ArgConvReverse_IndexClass.GetString(intCpnIndexType1);
		string cpnTerm1 = GetRefTerm1();
		int intCpnIndexDayCount1 = GetRefDayCount1();  
		string cpnIndexDayCount1 = ARM_ArgConvReverse_DayCount.GetString(intCpnIndexDayCount1);
		double refCoeff = GetRefCoeff1();

		//Number of columns
		size_t descSize = DEAL_NBCOLUMNS;
	
		size_t callEventSize = datesStructure.GetDateStrip(CALL_CRA_SCHED)->GetResetDates()->size();
		vector< string > rowDescVec(descSize);
		vector< ARM_GP_VALUE_TYPE > rowTypeVec(descSize, ARM_MISSING_TYPE); 
		
		InitPriceableColumns(rowDescVec,rowTypeVec);
			
		//Nb Non Call Fees
		int fwdStart;
		int j = 0;
		int feeSize = itsCallFees.size();
		double currFee = itsCallFees.GetDiscreteValues()->Elt(j);
		while ((currFee >= NON_CALL_FEE) && (j < feeSize-1))
		while ((j < feeSize) && (currFee >= NON_CALL_FEE))
		{
			j++;
			if (j < feeSize)
				currFee = itsCallFees.GetDiscreteValues()->Elt(j);
		}
		fwdStart = j;

		//NOTICE DATE
		double noticeDate;
		noticeDate = (*(datesStructure.GetDateStrip(CALL_CRA_SCHED)->GetResetDates()))[eventIdx];

		CC_Ostringstream noticeDateDesc;
		noticeDateDesc << CC_NS(std,fixed) << noticeDate;
		rowDescVec[ResetDate] = noticeDateDesc.str();
		rowTypeVec[ResetDate] = ARM_DATE_TYPE;

		//START DATE / CALL DATE
		double startDate =  (*(datesStructure.GetDateStrip(CALL_CRA_SCHED)->GetFlowStartDates()))[eventIdx] ;
		CC_Ostringstream startDateDesc;
		startDateDesc << CC_NS(std,fixed) << startDate;
		rowDescVec[StartDate] = startDateDesc.str();
		rowTypeVec[StartDate] = ARM_DATE_TYPE;
		
		//NEXT START DATE / NEXT CALL DATE
		double nextStartDate = (*(datesStructure.GetDateStrip(CALL_CRA_SCHED)->GetFlowEndDates()))[eventIdx] ;
		CC_Ostringstream nextStartDateDesc;
		nextStartDateDesc << CC_NS(std,fixed) << nextStartDate;
		rowDescVec[NextStartDate] = nextStartDateDesc.str();
		rowTypeVec[NextStartDate] = ARM_DATE_TYPE;

		//MATURITY DATE
		CC_Ostringstream maturityDateDesc;
		maturityDateDesc << CC_NS(std,fixed) << maturityDate;
		rowDescVec[MaturityDate] = maturityDateDesc.str();
		rowTypeVec[MaturityDate] = ARM_DATE_TYPE;

		//FEES
		double fees = const_cast< ARM_ReferenceValue& >(itsCallFees).Interpolate(startDate);
		CC_Ostringstream feesDesc;
		feesDesc << CC_NS(std,fixed) << fees;
		rowDescVec[Fees] = feesDesc.str();
		rowTypeVec[Fees] = ARM_DOUBLE;

		//FUNDING
		CC_Ostringstream fundingDesc;
		fundingDesc << "SWAP(" << ccy << ",";
		fundingDesc << "StartDate[i],";
		if(eventIdx == callEventSize-1)
			fundingDesc << "MaturityDate[i]" << ",";			// MaturityDate
		else
			fundingDesc << "NextStartDate[i],";					// FlowEndDate
		fundingDesc << "0," << payRec << ", ";
		fundingDesc << fundFreq << ",";
		fundingDesc << fundDayCount << ",";
		fundingDesc << fundFreq << ", ";
		fundingDesc << fundDayCount << ", ";
		fundingDesc << "FundSpreadCst,";
		fundingDesc << "NotionalCst" << ")";
		rowDescVec[Funding] = fundingDesc.str();
		rowTypeVec[Funding] = ARM_STRING;

		//CORRIDOR
		CC_Ostringstream corridorDesc;
		corridorDesc << "CORRIDOR(" << ccy << ",";			// Model
		corridorDesc << "StartDate[i]" << ",";				// StartDate
		if(eventIdx == callEventSize-1)
			corridorDesc << "MaturityDate[i]" << ",";		// MaturityDate
		else	
			corridorDesc << "NextStartDate[i]" << ",";		// FlowEndDate
		corridorDesc << payRec << ",";						// PayRec
		corridorDesc << boostedIndexType << ",";			// PayIndexType
		corridorDesc << "BoostedFixRateCst" << ",";			// Fix Value
		corridorDesc << 1				 << ",";			// PayIndexMult Value
		corridorDesc << "Z" << ",";							// Payment Freq
		corridorDesc << boostedVarTerm << ",";				// Pay Index Term
		corridorDesc << boostedResetTiming << ",";			// Pay Index Timing
		corridorDesc << boostedPayDayCount << ",";			// Pay Day Count
		corridorDesc << cpnResetGap << ",";					// Reset Gap
		corridorDesc << cpnIntRule << ",";					// Int Rule
		corridorDesc << "CpnSpreadCst" << ",";				// Spread
		corridorDesc << cpnResetFreq << ",";					// Fixing Freq
		corridorDesc << cpnResetTiming << ",";				// Fixing Timing
		corridorDesc << cpnIndexType1 << ",";				// Fixing Index Type 1
		corridorDesc << cpnTerm1 << ",";						// Fixing Index Term 1
		corridorDesc << refCoeff << ",";				// Coeff 1
		corridorDesc << "BarrierDownCst" << ",";				// Barrier Down
		corridorDesc << "BarrierUpCst" << ",";				// Barrier Up
		corridorDesc << "NotionalCst" << ")";				// Notional
		
		rowDescVec[Corridor] = corridorDesc.str();
		rowTypeVec[Corridor] = ARM_STRING;

		//EXOTIC SWAP
		CC_Ostringstream ExoticSwapDesc;
		double wCorridor = 0.0;
		double wFunding  = 0.0;
		if(ARM_CorridorLeg* corridorLeg = dynamic_cast<ARM_CorridorLeg*>(itsOptionPortfolio->GetPtf()->GetAsset(0)))
		{
			wCorridor = itsOptionPortfolio->GetPtf()->GetWeights()->Elt(0);
			wFunding  = itsOptionPortfolio->GetPtf()->GetWeights()->Elt(1);
		}
		else
		{
			wCorridor = itsOptionPortfolio->GetPtf()->GetWeights()->Elt(1);
			wFunding  = itsOptionPortfolio->GetPtf()->GetWeights()->Elt(0);
		}
	
		ExoticSwapDesc << wCorridor << "*Corridor[i]+" << wFunding << "*Funding[i]";
		rowDescVec[ExoticSwap] = ExoticSwapDesc.str();
		rowTypeVec[ExoticSwap] = ARM_STRING;

		//EXOTIC SWAP 1 
		CC_Ostringstream ExoticSwap1Desc;
		if(eventIdx == callEventSize-1)
		{
			ExoticSwap1Desc << "ExoticSwap[i]";
		}
		else
		{
			ExoticSwap1Desc << "ExoticSwap[i] + PV(" << "ExoticSwap1[i+1])";
		}
		rowDescVec[ExoticSwap1] = ExoticSwap1Desc.str();
		rowTypeVec[ExoticSwap1] = ARM_STRING;   
		
		//OPTION
		CC_Ostringstream optionDesc;
		if(eventIdx == callEventSize-1)
		{
			optionDesc << "MAX(" << "ExoticSwap1[i]" << "-" << "Fees[i]," << "0" << ")";
		}
		else
		{
			optionDesc << "MAX(" << "ExoticSwap1[i]" << "-" << "Fees[i]," << "PV(Option[i+1])" << ")";	
		}

		rowDescVec[Option] = optionDesc.str();
		rowTypeVec[Option] = ARM_STRING;

		//BERMUDA 1 
		CC_Ostringstream bermudaDesc;
		if (eventIdx == 0) 
		{
			if (itsOptionPortfolio->GetPorS() == K_PAY)
			{
				bermudaDesc << "-Option[i]";
			}
			else
			{
				bermudaDesc << "Option[i]";
			}
			rowDescVec[Bermuda] = bermudaDesc.str();
			rowTypeVec[Bermuda] = ARM_STRING;
		}
		else 
		{
			bermudaDesc << "0";	
			rowDescVec[Bermuda] = bermudaDesc.str();
			rowTypeVec[Bermuda] = ARM_DOUBLE;
		}

		//FUNDING LEG
		CC_Ostringstream fundingLegDesc;
		fundingLegDesc << "SWAP(" << ccy << ", StartDate[i], MaturityDate[i], 0," << payRec << ", ";
		fundingLegDesc << fundFreq << ", " << fundDayCount << ", " << fundFreq << ", " << fundDayCount << ", ";
		fundingLegDesc << "FundSpreadCst, NotionalCst" << ")";
		rowDescVec[FundingLeg] = fundingLegDesc.str();
		rowTypeVec[FundingLeg] = ARM_STRING;

		//CORRIDOR LEG
		CC_Ostringstream corridorLegDesc;
		corridorLegDesc << "CORRIDOR(" << ccy << ",";			// Model
		corridorLegDesc << "StartDate[i]" << ",";				// StartDate
		corridorLegDesc << "MaturityDate[i]" << ",";			// EndDate
		corridorLegDesc << payRec << ",";						// PayRec
		corridorLegDesc << boostedIndexType << ",";				// PayIndexType
		corridorLegDesc << "BoostedFixRateCst" << ",";			// Fix Value
		corridorLegDesc << 1					<< ",";			// PayIndexMult Value
		corridorLegDesc << cpnPayFreq << ",";					// Payment Freq
		corridorLegDesc << boostedVarTerm << ",";				// Pay Index Term
		corridorLegDesc << boostedResetTiming << ",";			// Pay Index Timing
		corridorLegDesc << boostedPayDayCount << ",";			// Pay Day Count
		corridorLegDesc << cpnResetGap << ",";					// Reset Gap
		corridorLegDesc << cpnIntRule << ",";					// Int Rule
		corridorLegDesc << "CpnSpreadCst" << ",";				// Spread
		corridorLegDesc <<  cpnResetFreq << ",";				// Fixing Freq
		corridorLegDesc << cpnResetTiming << ",";				// Fixing Timing
		corridorLegDesc << cpnIndexType1 << ",";				// Fixing Index Type 1
		corridorLegDesc << cpnTerm1 << ",";						// Fixing Index Term 1
		corridorLegDesc << refCoeff			<< ",";				// Coeff 1
		corridorLegDesc << "BarrierDownCst" << ",";				// Barrier Down
		corridorLegDesc << "BarrierUpCst" << ",";				// Barrier Up
		corridorLegDesc << "NotionalCst" << ")";				// Notional
		
		rowDescVec[CorridorLeg] = corridorLegDesc.str();
		rowTypeVec[CorridorLeg] = ARM_STRING;

		//EXOTIC SWAP 2 
		CC_Ostringstream ExoticSwap2Desc;
		ExoticSwap2Desc << wCorridor << "*CorridorLeg[i]+" << wFunding << "*FundingLeg[i]";
		rowDescVec[ExoticSwap2] = ExoticSwap2Desc.str();
		rowTypeVec[ExoticSwap2] = ARM_STRING;

		//OPTION 2  
		CC_Ostringstream option2Desc;
		if(eventIdx == callEventSize-1)
		{
			option2Desc << "MAX(ExoticSwap2[i]-Fees[i],0)";
		}
		else
		{
			option2Desc << "Exercise(0,ExoticSwap2[i]-Fees[i],Option2[i+1])";
		}
		rowDescVec[Option2] = option2Desc.str();
		rowTypeVec[Option2] = ARM_STRING;

		//BERMUDA 2 
		CC_Ostringstream bermuda2Desc;
		if(eventIdx == 0)
		{
			if (itsOptionPortfolio->GetPorS() == K_PAY)
			{
				bermuda2Desc << "-Option2[i]";
			}
			else
			{
				bermuda2Desc << "Option2[i]";
			}
			rowDescVec[Bermuda2] = bermuda2Desc.str();
			rowTypeVec[Bermuda2] = ARM_STRING;
		}
		else
		{
			bermuda2Desc << "0";
			rowDescVec[Bermuda2] = bermuda2Desc.str();
			rowTypeVec[Bermuda2] = ARM_DOUBLE;
		}

		//FWD CORRIDOR / FUNDING
		CC_Ostringstream fwdCorridorDesc;
		CC_Ostringstream fwdFundingDesc;
		if (eventIdx == fwdStart)
		{
			fwdCorridorDesc << "CorridorLeg[i]";
			rowDescVec[FwdCorridor] = fwdCorridorDesc.str();
			rowTypeVec[FwdCorridor] = ARM_STRING;
			
			fwdFundingDesc << "FundingLeg[i]";
			rowDescVec[FwdFunding] = fwdFundingDesc.str();
			rowTypeVec[FwdFunding] = ARM_STRING;
		}
		else
		{
			fwdCorridorDesc << "0";
			rowDescVec[FwdCorridor] = fwdCorridorDesc.str();
			rowTypeVec[FwdCorridor] = ARM_DOUBLE;

			CC_Ostringstream fwdFundingDesc;
			fwdFundingDesc << "0";
			rowDescVec[FwdFunding] = fwdFundingDesc.str();
			rowTypeVec[FwdFunding] = ARM_DOUBLE;
		}

		return ARM_RowInfo(rowDescVec,rowTypeVec);
	}
	catch (Exception& x)
	{
	    x.DebugPrint();

		throw x;
	}
	catch (...)
	{
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in CraCalculator::MiddleRows" );
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRACalculator
///	Routine: InitPriceableColumns
///	Returns: nothing
///	Action : Feed the MDM with input market objects
/////////////////////////////////////////////////////////////////
void ARM_CRACalculator::InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const
{
    string zeroValue("0");

	rowDescVec[Bermuda] = zeroValue;
	rowTypeVec[Bermuda] = ARM_DOUBLE;

	rowDescVec[Bermuda2] = zeroValue;
	rowTypeVec[Bermuda2] = ARM_DOUBLE;

	rowDescVec[Funding] = zeroValue;
	rowTypeVec[Funding] = ARM_DOUBLE;

	rowDescVec[Corridor] = zeroValue;
	rowTypeVec[Corridor] = ARM_DOUBLE;

	rowDescVec[FwdFunding] = zeroValue;
	rowTypeVec[FwdFunding] = ARM_DOUBLE;

	rowDescVec[FwdCorridor] = zeroValue;
	rowTypeVec[FwdCorridor] = ARM_DOUBLE;

	rowDescVec[ExoticSwap] = zeroValue;
	rowTypeVec[ExoticSwap] = ARM_DOUBLE;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRACalculator
///	Routine: ReplaceExerDates
///	Returns: void
///	Action : replace the reset dates of the date strip by
///          the exercise dates of the CRA and update fee table
///          to disable exercise at non relevant dates
/////////////////////////////////////////////////////////////////
std::vector<double>& ARM_CRACalculator::ComputeExerciseDatesAndSetFees() const
{
	ARM_Vector* StartDates = NULL;
	ARM_Vector* fees = NULL;
	std::vector<double>& NewResetDates = NULL;



	if (itsOptionPortfolio)
	{
		ARM_CorridorLeg* corridor = (ARM_CorridorLeg*)itsOptionPortfolio->GetCorridorLeg();
	
		StartDates = corridor->GetFlowStartDates();
		ARM_Vector* ResetDates = corridor->GetResetDates();
		
		int legSize = ResetDates->GetSize();

		ARM_Vector* NoticeDates = itsOptionPortfolio->GetExerciseStyle()->GetExerciseStartDates();
		ARM_Vector* ExpiryDates = itsOptionPortfolio->GetExerciseStyle()->GetExerciseEndDates();
		int noticeSize = NoticeDates->GetSize();

		NewResetDates = To_pARM_GP_Vector(ResetDates);

		fees = new ARM_Vector(legSize, NON_CALL_FEE);

		int lastIdx = 0;
		for (int i = 0; i < noticeSize; i++)
		{
			for (int j = 0; j < legSize; j++)
			{
				if (ARM_Date(ExpiryDates->Elt(i)) == ARM_Date(StartDates->Elt(j)))
				{
					NewResetDates->Elt(j) = NoticeDates->Elt(i);
					fees->Elt(j) = itsOptionPortfolio->GetStrike()->GetDiscreteValues()->Elt(i);
					break;
				}
			}
		
			if (j == legSize)
			{
				char msg[200];
				sprintf(msg, "ARM_CRACalculator::ComputeExerciseDatesAndSetFees : NoticeDates[%d] does not match any corridor start date", i);
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg);
			}
		}
	}
	else
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
						 "ARM_CRACalculator::ComputeExerciseDatesAndSetFees : The option portfolio and swaption should not be null.");
	}
	
	itsCallFees.SetDiscreteDates(StartDates);
	itsCallFees.SetDiscreteValues(fees);
	itsCallFees.SetCalcMethod(K_STEPUP_LEFT);

	return NewResetDates;
}




/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRACalculator
///	Routine: CustomDatesStructure
///	Returns: ARM_DateStripCombiner
///	Action : create the list of all event dates of the CRA. 
///			 The DateStripCombiner merges event dates of each
///          legs
///			customized Schedule Not implemented yet.
/////////////////////////////////////////////////////////////////
ARM_DateStripCombiner ARM_CRACalculator::CustomDatesStructure(const ARM_DateStripVector& dateStrips) const
{
	return ARM_DateStripCombiner();
}

/////////////////////////////////////////////////////////////////////////////
///	Class  : ARM_CRACalculator
///	Routine: DatesStructure
///	Returns: ARM_DateStripVector
///	Action : create the list of all event dates of the CRA. 
///			 The DateStripCombiner merges event dates of each
///          legs
/////////////////////////////////////////////////////////////////
ARM_DateStripCombiner ARM_CRACalculator::DatesStructure() const
{
	try
	{
		//General datas
		ARM_Currency ccy = GetCcy();
		
		ARM_DateStripVector SchedVect(1,NULL);
	
		// We take some call informations from the coupon leg.
		int dayCount = itsBoostedDayCount;
		int adjRule	= itsBoostedAdjRule;
		int rule = itsBoostedIntRule;
	
		// In this case the exercise dates are computed from the 
		// option portfolio data
		if (itsIsPortfolioNoticeDays)
		{
			if (itsOptionPortfolio)
			{
				ARM_CorridorLeg* corridor = (ARM_CorridorLeg*)itsOptionPortfolio->GetCorridorLeg();
	
				std::vector<double>& FlowStartDates		= new std::vector<double>();
				std::vector<double>& FlowEndDates			= new std::vector<double>();
				std::vector<double>& FwdRateStartDates	= new std::vector<double>();
				std::vector<double>& FwdRateEndDates		= new std::vector<double>();
				std::vector<double>& PaymentDates			= new std::vector<double>();
				std::vector<double>& InterestTerms		= new std::vector<double>();
				std::vector<double>& ResetDates			= ComputeExerciseDatesAndSetFees();
				// NB : from CCSO, we may have less notice dates than corridor's periods (notice freq < payment freq)
				// in that case, we skip the periods that do not have a notice (no exercise => PV = 0)
				// these periods were set with a reset date = -1
				std::vector<double>& FlowResetDates		= new std::vector<double>();

				int size = corridor->GetFlowStartDates()->size();
				
				std::vector<double>& InterestDays = new std::vector<double>();
				for (int i=0; i<size; i++)
				{
					if (ResetDates->Elt(i) != -1.0)
					{
						FlowStartDates->push_back(corridor->GetFlowStartDates()->Elt(i));
						FlowEndDates->push_back(corridor->GetFlowEndDates()->Elt(i));
						FwdRateStartDates->push_back(corridor->GetPaymentStartDates()->Elt(i));
						FwdRateEndDates->push_back(corridor->GetPaymentEndDates()->Elt(i));
						PaymentDates->push_back(corridor->GetPaymentDates()->Elt(i));
						InterestTerms->push_back(corridor->GetPaymentInterestterms()->Elt(i));
						FlowResetDates->push_back(ResetDates->Elt(i));

						InterestDays->push_back(DaysBetweenDates(dayCount,
																 FlowStartDates->Elt(i), 
																 FlowEndDates->Elt(i)));
					}
				}
				
				ARM_DateStrip CallSched(FlowStartDates, FlowEndDates, FwdRateStartDates, FwdRateEndDates,
										FlowResetDates, PaymentDates, InterestDays, InterestTerms);

				SchedVect[CALL_CRA_SCHED] = &CallSched;

				delete FlowStartDates;
				FlowStartDates = NULL;

				delete FlowEndDates;
				FlowEndDates = NULL;
				
                delete ResetDates;
				ResetDates = NULL;
				
                delete FwdRateStartDates;
				FwdRateStartDates = NULL;
				
                delete FwdRateEndDates;
				FwdRateEndDates = NULL;
				
                delete PaymentDates;
				PaymentDates = NULL;
				
                delete InterestDays;
				InterestDays = NULL;
				
                delete InterestTerms;
				InterestTerms = NULL;
				
				delete ResetDates;
				ResetDates = NULL;

				delete FlowResetDates;
				FlowResetDates = NULL;

				itsExerDateStrip = ARM_DateStripPtr(static_cast<ARM_DateStrip*>(CallSched.Clone()));

				return SchedVect;
			}
			else if (itsSwaption)
			{
				ARM_SwapLeg* fundingLeg			= (ARM_SwapLeg*) itsSwaption->GetFloatLeg();
				ARM_Vector* vFlowStartDates		= fundingLeg->GetFlowStartDates();
				ARM_Vector* vFlowEndDates		= fundingLeg->GetFlowEndDates();
				ARM_Vector* vResetDates			= fundingLeg->GetResetDates();
				ARM_Vector* vFwdRateStartDates	= fundingLeg->GetFwdRateStartDates();
				ARM_Vector* vFwdRateEndDates	= fundingLeg->GetFwdRateEndDates();
				ARM_Vector* vPaymentDates		= fundingLeg->GetPaymentDates();
				ARM_Vector* vInterestTerms		= fundingLeg->GetInterestTerms();

				ARM_Vector* NoticeDates = itsSwaption->GetExerciseStyle()->GetExerciseStartDates();
				ARM_Vector* ExpiryDates = itsSwaption->GetExerciseStyle()->GetExerciseEndDates();

				int noticeSize = NoticeDates->GetSize();
				int legSize    = vResetDates->GetSize();

				std::vector<double>& FlowStartDates		= new std::vector<double>(noticeSize);
				std::vector<double>& FlowEndDates			= new std::vector<double>(noticeSize);
				std::vector<double>& ResetDates			= new std::vector<double>(noticeSize);
				std::vector<double>& FwdRateStartDates	= new std::vector<double>(noticeSize);
				std::vector<double>& FwdRateEndDates	    = new std::vector<double>(noticeSize);
				std::vector<double>& PaymentDates			= new std::vector<double>(noticeSize);
				std::vector<double>& InterestTerms		= new std::vector<double>(noticeSize);
				std::vector<double>& InterestDays			= new std::vector<double>(noticeSize);

				ARM_Vector* fees = NULL;
				ARM_Vector* StartDatesfees = new ARM_Vector(noticeSize);

				for (int i = 0; i < noticeSize; i++)
				{
					for (int j = 0; j < legSize; j++)
					{
						if (ARM_Date(ExpiryDates->Elt(i)) == ARM_Date(vFlowStartDates->Elt(j)))
						{
							ResetDates->Elt(i)		  = NoticeDates->Elt(i);
							FlowStartDates->Elt(i)    = vFlowStartDates->Elt(j);
							StartDatesfees->Elt(i)    = FlowStartDates->Elt(i);
							FlowEndDates->Elt(i)	  = vFlowEndDates->Elt(j);
							FwdRateStartDates->Elt(i) = vFwdRateStartDates->Elt(j);
							FwdRateEndDates->Elt(i)   = vFwdRateEndDates->Elt(j);
							PaymentDates->Elt(i)      = vPaymentDates->Elt(j);
							InterestTerms->Elt(i)     = vInterestTerms->Elt(j);
							break;
						}
					}
					if (j == legSize)
					{
						char msg[200];
						sprintf(msg, "ARM_CRACalculator::ComputeExerciseDatesAndSetFees : NoticeDates[%d] does not match any corridor start date", i);
						throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg);
					}
				}
				
				for (i=0; i<noticeSize; i++)
				{
					InterestDays->Elt(i) = DaysBetweenDates(dayCount,
															FlowStartDates->Elt(i), 
															FlowEndDates->Elt(i));
				}

				ARM_DateStrip CallSched(FlowStartDates, FlowEndDates, FwdRateStartDates, FwdRateEndDates,
										ResetDates, PaymentDates, InterestDays, InterestTerms);

				SchedVect[CALL_CRA_SCHED] = &CallSched;

				// process fees
				if (itsSwaption->GetFee() != NULL)
				{
					//Adjust summit fees
					int feeSize = itsSwaption->GetFee()->size();
					ARM_ReferenceValue newFees;
					ARM_Vector* vecNoticeDates	= itsSwaption->GetFee()->GetDiscreteDates();
					ARM_Vector* vecFees			= itsSwaption->GetFee()->GetDiscreteValues();
					ARM_Vector* vecNewFees		= NULL; 
					double sumFee				= vecFees->Elt(0);
					double firstFee				= 0.0;
					double nextFee				= 0.0;
					bool isZero					= false;
					bool isOne					= false;

					for (i = 1; i < feeSize; i++)
					{
						firstFee = vecFees->Elt(i-1);
						nextFee  = vecFees->Elt(i);
						sumFee += nextFee;

						if ((firstFee == 1) && (firstFee == nextFee))
						{
							isOne = true;
						}
						else
						{
							isOne = false;
						}
					}
					
					//Keep the current fees
					if (sumFee == 0.0)
					{
						isZero = true;
						itsCallFees = ARM_ReferenceValue((ARM_Vector*)vecNoticeDates->Clone(), (ARM_Vector*) vecFees->Clone());
					}
					//Modify the fees from ONE to ZERO
					else if (isOne == true)
					{
						for (i = 0; i < feeSize; i++)
						{
							vecNewFees = new ARM_Vector(feeSize, 0.0); 
							itsCallFees = ARM_ReferenceValue((ARM_Vector*)vecNoticeDates->Clone(), (ARM_Vector*) vecNewFees->Clone());
						}
					}
					//Custom Fee
					else
					{
						vecNewFees = new ARM_Vector(feeSize, 0.0); 
						double currFee = 0.0;
						double currNoticeDate = 0.0;
						double currNominal = 0.0;
						for (i = 0; i < feeSize; i++)
						{
							currFee = vecFees->Elt(i);
							currNoticeDate = vecNoticeDates->Elt(i);
							if ( currFee == 1.0)
							{
								vecNewFees->Elt(i) = 0.0;
							}
							else if (currFee == 0.0)
							{
								vecNewFees->Elt(i) = 0.0;
							}
							else
							{
								if (itsSwaption->GetAmount()->size() == 1)
									currNominal = itsSwaption->Get1stLeg()->GetAmount()->GetDiscreteValues()->Elt(0);
								else
									currNominal = itsSwaption->Get1stLeg()->GetAmount()->Interpolate(currNoticeDate);

								if (itsSwaption->GetFixedLeg()->GetRcvOrPay() == K_PAY)
								{
									if (currFee < 1.0)
										vecNewFees->Elt(i) = (1.0 - currFee)*currNominal;				
									else //(currFee > 1) 
										vecNewFees->Elt(i) = -(currFee - 1.0)*currNominal;
								}
								else
								{
									if (currFee < 1.0)
										vecNewFees->Elt(i) = -(1.0 - currFee)*currNominal;				
									else //(currFee > 1) 
										vecNewFees->Elt(i) = (currFee - 1.0)*currNominal;
								}
							}
						}

						itsCallFees = ARM_ReferenceValue((ARM_Vector*)vecNoticeDates->Clone(), (ARM_Vector*) vecNewFees->Clone());
					}
					
					if (vecNewFees)
					{
						delete vecNewFees;
						vecNewFees = NULL;
					}
				}
				else
				{
					ARM_Vector* fees = new ARM_Vector(noticeSize,0.0);
					itsCallFees.SetDiscreteDates(StartDatesfees);
					itsCallFees.SetDiscreteValues(fees);
				}
				itsCallFees.SetCalcMethod(K_STEPUP_LEFT);

				delete FlowStartDates;
				FlowStartDates = NULL;

				delete FlowEndDates;
				FlowEndDates = NULL;
				
                delete ResetDates;
				ResetDates = NULL;
				
                delete FwdRateStartDates;
				FwdRateStartDates = NULL;
				
                delete FwdRateEndDates;
				FwdRateEndDates = NULL;
				
                delete PaymentDates;
				PaymentDates = NULL;
				
                delete InterestDays;
				InterestDays = NULL;
				
                delete InterestTerms;
				InterestTerms = NULL;

				delete StartDatesfees;
				StartDatesfees = NULL;
				
				itsExerDateStrip = ARM_DateStripPtr(static_cast<ARM_DateStrip*>(CallSched.Clone()));

				return SchedVect;
			}
			else
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
								 "ARM_CRACalculator::DatesStructure : The option portfolio and swaption should not be null.");
			}
		}
		// The dates cames from the datestrip
		else
		{
			//Call schedule
			int callFreq				= GetCallFreq();
			int callNotice				= GetCallNotice();
			const char* callCal			= GetCallCal();
			
			ARM_ReferenceValue callFees = itsCallFees;
			double first				= ((*callFees.GetDiscreteDates()))[0];
			ARM_Date firstCall			= ARM_Date(first);
			int size					= callFees.size();
			double last					= ((*callFees.GetDiscreteDates()))[size-1];

			double startDate, endDate;
		
			startDate = itsStartDate.GetJulian();
			endDate = itsEndDate.GetJulian();
	
			ARM_DateStrip CallSched(startDate, endDate, callFreq, dayCount,
									callCal, adjRule, rule, K_SHORTSTART, callNotice,
									callFreq, GETDEFAULTVALUE, callCal, K_ADVANCE, K_ARREARS, true);

			SchedVect[CALL_CRA_SCHED]    = &CallSched;

			itsExerDateStrip = ARM_DateStripPtr(static_cast<ARM_DateStrip*>(CallSched.Clone()));

			return SchedVect;
		}
	}

	catch (Exception& x)
	{
	    x.DebugPrint();

		throw x;
	}

	catch (...)
	{
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in CraCalculator::DatesStructure" );
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRACalculator
///	Routine: GenerateOptionPortfolio
///	Returns: ARM_OptionPortfolio*
///	Action : creates the instrument (swap + corridor)
/////////////////////////////////////////////////////////////////
ARM_OptionPortfolio* ARM_CRACalculator::GenerateOptionPortfolio()
{
	try
	{
		ARM_OptionPortfolio* optionPort = NULL;
		ARM_CorridorLeg* corridorLeg = NULL;
		ARM_SwapLeg* fundingLeg = NULL;
		
		//Useful datas
		ARM_Date startDate = GetStartDate();
		ARM_Date endDate = GetEndDate();
		int payRec = GetPayRec();
		int fundPayRec = ((payRec == K_RCV) ? K_PAY : K_RCV);
		int refResetFreq = GetCpnResetFreq();
		int boostedResetTiming  = GetBoostedResetTiming();
		int cpnPayFreq= GetCpnPayFreq();
		
		ARM_IRIndex* paymentIndex = NULL; 
		ARM_ReferenceValue cpnSpread = GetCpnSpread();
		if (itsBoostedIndexType == Fixed)
		{
			paymentIndex = new ARM_IRIndex(GetCcy().GetCcyName(), itsBoostedDayCount);
			cpnSpread = ARM_ReferenceValue(itsBoostedFixRate);
		}
		else
		{
			string liborTypeName(string(GetCcy().GetCcyName()) == string("EUR") ? "EURIBOR" : "LIBOR");
			liborTypeName += itsBoostedVarTerm;
			ARM_INDEX_TYPE payIndexType = static_cast< ARM_INDEX_TYPE > (ARM_ArgConv_IndexType.GetNumber(liborTypeName));

			paymentIndex = new ARM_IRIndex( payIndexType, 
											K_DEF_FREQ, 
											K_DEF_FREQ, 
											&ARM_Currency(GetCcy().GetCcyName()), 
											itsBoostedDayCount);
			if (cpnPayFreq != K_DEF_FREQ)
				paymentIndex->SetPayFrequency(cpnPayFreq);
			paymentIndex->SetResetTiming(boostedResetTiming);
			paymentIndex->SetPayTiming(K_ARREARS);
		}
		
		cpnSpread *= 100.0;

		ARM_IRIndex* refIndex =NULL;
		string liborTypeName(string(GetCcy().GetCcyName()) == string("EUR") ? "EURIBOR" : "LIBOR");
		liborTypeName += itsRefTerm1;
		ARM_INDEX_TYPE refIndexType = static_cast< ARM_INDEX_TYPE > (ARM_ArgConv_IndexType.GetNumber(liborTypeName));

		refIndex = new ARM_IRIndex( refIndexType, 
									K_DEF_FREQ, 
									K_DEF_FREQ, 
									&ARM_Currency(GetCcy().GetCcyName()), 
									itsRefDayCount1);

		if (refResetFreq != K_DEF_FREQ)
			refIndex->SetResetFrequency(refResetFreq);
		
		int cpnResetTiming = GetCpnResetTiming();
		refIndex->SetResetTiming(cpnResetTiming);

		ARM_ReferenceValue cpnBarrierDown = GetCpnBarrierDown();
		cpnBarrierDown *= 100.0;
		ARM_ReferenceValue cpnBarrierUp = GetCpnBarrierUp();
		cpnBarrierUp *= 100.0;

   		ARM_Currency ccy = GetCcy();
		ARM_ReferenceValue fundSpread = GetFundSpread();
		fundSpread *= 100.0;

		ARM_ReferenceValue notional = GetNotional();
		notional.SetCalcMethod(K_STEPUP_RIGHT);

		int fundFreq = GetFundFreq();
		string fundTerm = ARM_ArgConvReverse_MatFrequency.GetString(fundFreq);
		
		liborTypeName = string(GetCcy().GetCcyName()) == string("EUR") ? "EURIBOR" : "LIBOR";
		if(fundTerm == "12M")
			fundTerm = "1Y"; // because we cant convert EURIBOR1Y and not EURIBOR12M
		liborTypeName += fundTerm;
		ARM_INDEX_TYPE fundIndex = static_cast< ARM_INDEX_TYPE > (ARM_ArgConv_IndexType.GetNumber(liborTypeName));
		
		//Portfolio
		corridorLeg = new ARM_CorridorLeg(startDate, 
										  endDate,
										  payRec, 
										  paymentIndex,
										  cpnPayFreq,
										  &cpnSpread,
										  refIndex, 
										  refResetFreq,
										  boostedResetTiming,	//PayIndex Reset Timing: case of libor pay index.
										  cpnResetTiming,		//Ref Index Reset Timing
										  K_SHORTSTART,			//TMP
										  &cpnBarrierDown, 
										  K_STD,
										  &cpnBarrierUp, 
										  K_STD,
										  &ccy,
										  K_DEF_FREQ, 
										  K_LINEAR,
										  K_DIGITALE,
										  0, //TMP: forced CdecompPricingFlag,
										  (char*)(GetCpnResetCal()).c_str(),
										  (char*)(GetCpnPayCal()).c_str());
		
		delete refIndex;
		refIndex = NULL;

		delete paymentIndex;
		paymentIndex = NULL;

		corridorLeg->SetAmount(&notional);

		fundingLeg =  new ARM_SwapLeg(startDate, 
									  endDate, 
									  (ARM_INDEX_TYPE) fundIndex, 
									  fundPayRec, 
									  0.0, 
									  GetFundFreq(), 
									  GetFundFreq(), 
									  K_ADVANCE, 
									  K_ARREARS,
									  &ccy,
									  K_ADJUSTED,
									  -ccy.GetSpotDays(),
									  (char*)(GetCpnResetCal()).c_str(),
									  (char*)(GetCpnPayCal()).c_str(),
									  1,
									  0,
									  K_SHORTSTART,
									  NULL,
									  1);

		fundingLeg->SetAmount(&notional);

		ARM_ReferenceValue theSpread = GetFundSpread();
        theSpread *= 100.0;
		theSpread.SetCalcMethod(K_STEPUP_LEFT);
		fundingLeg->SetVariableSpread(&theSpread);

		int nbsec=2;
		double* weights = new double[nbsec];
		double* mktprices = new double[nbsec];

		weights[0] = 1.0;
		weights[1] = 1.0;
		mktprices[0] = 1.0;
		mktprices[1] = 1.0;

		ARM_StdPortfolio* port = new ARM_StdPortfolio(nbsec, weights, mktprices);
		port->SetAsset(corridorLeg, 0);
		port->SetAsset(fundingLeg, 1);

		//Exercise Datas
		ARM_Vector tempExer(*GetCallFees().GetDiscreteDates());
		ARM_ReferenceValue genCallFees = itsCallFees;

		int sz = tempExer.GetSize();
		ARM_ExerciseStyle exerciseStyle(&tempExer);
		for (int i = 0; i < sz; i++)
		{
			ARM_Date tmpDate(exerciseStyle.GetExerciseStartDates()->Elt(i));
			tmpDate.PreviousBusinessDay(-itsCallNotice, ccy.GetCcyName());
			exerciseStyle.GetExerciseStartDates()->Elt(i) = tmpDate.GetJulian();
			exerciseStyle.GetExerciseEndDates()->Elt(i) = tmpDate.GetJulian();
			genCallFees.GetDiscreteDates()->Elt(i) = tmpDate.GetJulian();
		}

		//Option Portfolio
		optionPort = new ARM_OptionPortfolio(port, &exerciseStyle, &genCallFees, K_CALL);
		optionPort->SetPorS(payRec);

		//Free memory
		delete [] weights;
		weights = NULL;
		delete [] mktprices;
		mktprices = NULL;

		return optionPort;
	}
	catch (Exception& x)
	{
	    x.DebugPrint();

		throw x;
	}
	catch (...)
	{
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in CraCalculator::GenerateOptionPortfolio" );
	}
}

////////////////////////////////////////////////////////////////////////////
///	Class  : ARM_CRACalculator	
///	Routine: GenerateProductDescription	
///	Returns: void						
///	Action : creates the product description from summit option portfolio
////////////////////////////////////////////////////////////////////////////
void ARM_CRACalculator::GenerateProductDescription(ARM_OptionPortfolio* optPortfolio)
{
	try
	{
		int i=0;
		int size = 0;
		SetOptionPortfolio(optPortfolio);

		ARM_CorridorLeg* corridorLeg = optPortfolio->GetCorridorLeg();
		ARM_SwapLeg* fundingLeg = optPortfolio->GetLiborLeg();
		
		ARM_Vector* cpnStartDates = corridorLeg->GetFlowStartDates(); 
		ARM_Vector* cpnEndDates = corridorLeg->GetFlowEndDates(); 
		ARM_Vector* fundStartDates = fundingLeg->GetFlowStartDates(); 

		int cpnSize = cpnStartDates->size();

		//StartDate
		ARM_Date startDate = fundingLeg->GetStartDateNA();
		SetStartDate(startDate);
		
		//EndDate
		ARM_Date endDate = fundingLeg->GetEndDateNA();
		SetEndDate(endDate);
		
		//PayRec
		int payRec = corridorLeg->GetRcvOrPay();
		SetPayRec(payRec);
		
		//CouponNotional
		ARM_Vector* vecNotional = (ARM_Vector*)(corridorLeg->GetAmount()->GetDiscreteValues())->Clone();
		size = vecNotional->size();
		ARM_ReferenceValue notional;
		if (size == 1)
		{
			ARM_Vector* vecDates = new ARM_Vector(1, 0.0);
			vecDates->Elt(0) = startDate.GetJulian();
			notional = ARM_ReferenceValue(vecDates, vecNotional);
		}
		else
		{
			notional = *((ARM_ReferenceValue*)(corridorLeg->GetAmount())->Clone());
		}
		SetNotional(notional);
		
		//FundFreq
		int fundFreq = fundingLeg->GetPaymentFreq();
		SetFundFreq(fundFreq);
		
		//FundSpread
		if (fundingLeg->GetVariableSpread())
		{
			ARM_ReferenceValue* tempVarFundSpread = (ARM_ReferenceValue *) fundingLeg->GetVariableSpread()->Clone();
			*(tempVarFundSpread->GetDiscreteValues()) /= 100.0; 
			SetFundSpread(*tempVarFundSpread);
			delete tempVarFundSpread;
		}
		else
		{
			double spread = (fundingLeg->GetSpread())/100;
			ARM_Vector* vecDates = new ARM_Vector(1, startDate.GetJulian());
			ARM_Vector* vecFundSpread = new ARM_Vector(1, spread);
			ARM_ReferenceValue fundSpread = ARM_ReferenceValue(vecDates, vecFundSpread);
			SetFundSpread(fundSpread);
		}	
		
		//FundDayCount
		int fundDayCount = fundingLeg->GetIndexDayCount();
		SetFundDayCount(fundDayCount);
		
		//FundAdjRule
		itsFundAdjRule = fundingLeg->GetIRIndex()->GetFwdRule();
		
		//FundIntRule
		itsFundIntRule = fundingLeg->GetIRIndex()->GetIntRule();
		
		//FundStubRule
		itsFundStubRule = fundingLeg->GetStubMeth();
		
		//CouponIndexes
		ARM_IRIndex* cpnPayIndex = corridorLeg->GetPaymentIndex();
		ARM_IRIndex* cpnRefIndex = corridorLeg->GetRefIndex();
		
		//CpnPayFreq
		int cpnPayFreq = cpnPayIndex->GetPayFrequency();
		SetCpnPayFreq(cpnPayFreq);
		
		//CpnResetCal
		char* resetCal = corridorLeg->GetResetCalName();
		string cpnResetCal(resetCal); 
		SetCpnResetCal(cpnResetCal); 
		
		//CpnPayCal
		char* payCal = corridorLeg->GetPayCalName();
		string cpnPayCal(payCal);
		SetCpnPayCal(cpnPayCal); 
		
		//CpnResetGap
		int cpnResetGap = 0;
		SetCpnResetGap(cpnResetGap);
		
		//FundAdjRule
		itsCpnAdjRule = cpnRefIndex->GetFwdRule();
		
		//FundIntRule
		itsCpnIntRule = cpnRefIndex->GetIntRule();
		
		//FundStubRule
		itsCpnStubRule = corridorLeg->GetStubMeth();
		
		//BoostedIndexType & BoostedFixRate & CpnSpread
		int boostedIndexType;
		if (cpnPayIndex->IsFixedIndex())
		{
			boostedIndexType = Fixed;
			int period = 12/GetCcy().GetFixedPayFreq();
			if ( period == 12 )
			{
				itsBoostedVarTerm = string("1Y");
			}
			else
			{
				char str[2];
				sprintf(str, "%d", period);
				itsBoostedVarTerm = string(str) + string("M");
			}
			
			if (corridorLeg->GetSpreads())
			{
				ARM_ReferenceValue* boostedStepUpFixRate = (ARM_ReferenceValue *) corridorLeg->GetSpreads()->Clone();
				*(boostedStepUpFixRate->GetDiscreteValues()) /= 100.0;
				SetBoostedFixRate(*boostedStepUpFixRate);
				delete boostedStepUpFixRate;
			}
			else 
			{	
				double fixRate = (corridorLeg->GetSpread())/100;
				ARM_Vector* vecDates = new ARM_Vector(1, startDate.GetJulian());
				ARM_Vector* vecBoostedFixRate = new ARM_Vector(1, fixRate);
				ARM_ReferenceValue boostedFixRate = ARM_ReferenceValue(vecDates, vecBoostedFixRate);
				SetBoostedFixRate(boostedFixRate);
			}
		
			ARM_ReferenceValue* cpnSpread = new ARM_ReferenceValue(0.0);
			SetCpnSpread(*cpnSpread);
			delete cpnSpread;
		}
		else
		{
			boostedIndexType = Libor;

			if (corridorLeg->GetSpreads())
			{
				ARM_ReferenceValue* cpnStepUpSpread = (ARM_ReferenceValue *) corridorLeg->GetSpreads()->Clone();
				*(cpnStepUpSpread->GetDiscreteValues()) /= 100.0;
				SetCpnSpread(*cpnStepUpSpread);
				delete cpnStepUpSpread;
			}
			else
			{
				double spread = (corridorLeg->GetSpread())/100;
				ARM_Vector* vecDates = new ARM_Vector(1, startDate.GetJulian());
				ARM_Vector* vecCpnSpread = new ARM_Vector(1, spread);
				ARM_ReferenceValue cpnSpread = ARM_ReferenceValue(vecDates, vecCpnSpread);
				SetCpnSpread(cpnSpread);
			}

			ARM_ReferenceValue* boostedFixedRate = (ARM_ReferenceValue*)itsCpnSpread.Clone();
			SetBoostedFixRate(*boostedFixedRate);
			delete boostedFixedRate;
		}
		SetBoostedIndexType(boostedIndexType);
		
		//BoostedVarTerm
		string boostedVarTerm;
		int intBoostedVarTerm = cpnPayIndex->GetTerm();
		if (intBoostedVarTerm == K_ANNUAL)
			boostedVarTerm = "1Y";
		else if (intBoostedVarTerm == K_SEMIANNUAL)
			boostedVarTerm = "6M";
		else if (intBoostedVarTerm == K_QUARTERLY)
			boostedVarTerm = "3M";
		else if (intBoostedVarTerm == K_MONTHLY)
			boostedVarTerm = "1M";
		SetBoostedVarTerm(boostedVarTerm);
		
		//BoostedResetGap
		int boostedResetGap = cpnPayIndex->GetResetGap();
		SetBoostedResetGap(boostedResetGap);
		
		//BoostedResetTiming
		int boostedResetTiming = cpnPayIndex->GetResetTiming();
		SetBoostedResetTiming(boostedResetTiming);

		//BoostedDayCount
		int boostedDayCount = cpnPayIndex->GetDayCount(); 
		SetBoostedDayCount(boostedDayCount);
		
		//BoostedAdjRule
		int boostedAdjRule = cpnPayIndex->GetFwdRule();
		SetBoostedAdjRule(boostedAdjRule);
		
		//BoostedIntRule
		int boostedIntRule = cpnPayIndex->GetIntRule();
		SetBoostedIntRule(boostedIntRule);
		
		//CpnBarrierDown
		ARM_ReferenceValue* cpnBarrierDown = (ARM_ReferenceValue *) corridorLeg->GetDownBarriers()->Clone();
		*(cpnBarrierDown->GetDiscreteValues()) /= 100.0;
		SetCpnBarrierDown(*cpnBarrierDown);
		delete cpnBarrierDown;
		
		//CpnBarrierUp
		ARM_ReferenceValue* cpnBarrierUp = (ARM_ReferenceValue *) corridorLeg->GetUpBarriers()->Clone();
		*(cpnBarrierUp->GetDiscreteValues()) /= 100.0;
		SetCpnBarrierUp(*cpnBarrierUp);
		delete cpnBarrierUp;
				
		//CpnResetFreq
		int cpnResetFreq = cpnRefIndex->GetResetFrequency(); //K_MONTHLY;
		SetCpnResetFreq(cpnResetFreq);
		
		//CpnResetTiming
		int cpnResetTiming = cpnRefIndex->GetResetTiming();
		SetCpnResetTiming(cpnResetTiming);
		
		//NB: One of the two references indexes is not used.
		//RefIndexType
		int refIndexType;
		if (cpnRefIndex->IsLiborIndex())
			refIndexType = Libor;
		else
		{
			if (cpnRefIndex->IsCMSIndex())
				refIndexType = Cms;
		}
		SetRefIndexType1(refIndexType);
		
		//RefTerm
		string refTerm; 
		int intTerm = cpnRefIndex->GetTerm();
		if (intTerm == K_ANNUAL)
			refTerm = "1Y";
		else if (intTerm == K_SEMIANNUAL)
			refTerm = "6M";
		else if (intTerm == K_QUARTERLY)
			refTerm = "3M";
		else if (intTerm == K_MONTHLY)
			refTerm = "1M";
		SetRefTerm1(refTerm);
		
		//RefDayCount
		int refDayCount = cpnRefIndex->GetDayCount();
		SetRefDayCount1(refDayCount);
		
		SetRefCoeff1(1);
	
		//CallCal
		string callCal = cpnResetCal;
		SetCallCal(callCal);
		char* tempCallCal = (char *) callCal.c_str();
		
		//CallFees
		//Call notice : no need when we are in Summit case
		SetCallNotice(9999);

		//CallFreq
		int callFreq = cpnPayFreq; 
		SetCallFreq(callFreq);
	}
	catch (Exception& x)
	{
	    x.DebugPrint();

		throw x;
	}
	catch (...)
	{
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in CraCalculator::GenerateProductDescription" );
	}
}

////////////////////////////////////////////////////////////////////////////
///	Class  : ARM_CRACalculator	
///	Routine: GenerateProductDescription	
///	Returns: void						
///	Action : creates the product description from summit swaption
////////////////////////////////////////////////////////////////////////////
void ARM_CRACalculator::GenerateProductDescription(ARM_Swaption* swaption)
{
	try
	{
		int i=0;
		int size = 0;
		SetSwaption(swaption);

		ARM_SwapLeg* fundingLeg = itsSwaption->GetFloatLeg();
		ARM_SwapLeg* fixedLeg = itsSwaption->GetFixedLeg();
		
		ARM_Vector* fundStartDates = fundingLeg->GetFlowStartDates(); 
		int cpnSize = fundStartDates->size();

		//StartDate
		ARM_Date startDate = fundingLeg->GetStartDateNA();
		SetStartDate(startDate);
		
		//EndDate
		ARM_Date endDate = fundingLeg->GetEndDateNA();
		SetEndDate(endDate);
		
		//PayRec
		int payRec = itsSwaption->GetOptionType();
		SetPayRec(payRec);
		
		//CouponNotional
		size = fundingLeg->GetAmount()->size();
		ARM_ReferenceValue notional;
		if (size == 1)
		{
			ARM_Vector* vecValues = (ARM_Vector*)(fundingLeg->GetAmount()->GetDiscreteValues())->Clone();
			ARM_Vector* vecDates = new ARM_Vector(1, startDate.GetJulian());
			notional = ARM_ReferenceValue(vecDates, vecValues);
		}
		else
		{
			notional = *((ARM_ReferenceValue*)(fundingLeg->GetAmount())->Clone());
		}
		SetNotional(notional);
		
		//FundFreq
		int fundFreq = fundingLeg->GetPaymentFreq();
		SetFundFreq(fundFreq);
		
		//FundSpread
		if (fundingLeg->GetVariableSpread())
		{
			ARM_ReferenceValue* tempVarFundSpread = (ARM_ReferenceValue*) fundingLeg->GetVariableSpread()->Clone();
			*(tempVarFundSpread->GetDiscreteValues()) /= 100.0; 
			SetFundSpread(*tempVarFundSpread);
			delete tempVarFundSpread;
		}
		else
		{
			double spread = (fundingLeg->GetSpread())/100;
			ARM_Vector* vecDates = new ARM_Vector(1, startDate.GetJulian());
			ARM_Vector* vecFundSpread = new ARM_Vector(1, spread);
			ARM_ReferenceValue fundSpread = ARM_ReferenceValue(vecDates, vecFundSpread);
			SetFundSpread(fundSpread);
		}	
		
		//FundDayCount
		int fundDayCount = fundingLeg->GetIndexDayCount();
		SetFundDayCount(fundDayCount);
		
		//FundAdjRule
		itsFundAdjRule = fundingLeg->GetIRIndex()->GetFwdRule();
		
		//FundIntRule
		itsFundIntRule = fundingLeg->GetIRIndex()->GetIntRule();
		
		//FundAStubRule
		itsFundStubRule = fundingLeg->GetStubMeth();
		
		//CouponIndexes
		ARM_IRIndex* cpnPayIndex = fixedLeg->GetIRIndex();
		ARM_IRIndex* cpnRefIndex = fundingLeg->GetIRIndex();
		
		//CpnPayFreq
		int cpnPayFreq = fundingLeg->GetPaymentFreq();
		SetCpnPayFreq(cpnPayFreq);
		
		//CpnResetCal
		char* resetCal = fundingLeg->GetResetCalName();
		string cpnResetCal(resetCal); 
		SetCpnResetCal(cpnResetCal); 
		
		//CpnPayCal
		char* payCal = fundingLeg->GetPayCalName();
		string cpnPayCal(payCal);
		SetCpnPayCal(cpnPayCal); 
		
		//CpnResetGap
		int cpnResetGap = 0;
		SetCpnResetGap(cpnResetGap);
		
		//CpnAdjRule
		itsCpnAdjRule = fundingLeg->GetIRIndex()->GetFwdRule();
		
		//CpnIntRule
		itsCpnIntRule = fundingLeg->GetIRIndex()->GetIntRule();
		
		//CpnAStubRule
		itsCpnStubRule = fundingLeg->GetStubMeth();

		//BoostedIndexType & BoostedFixRate & CpnSpread
		int boostedIndexType = Fixed;
		int period = 12/GetCcy().GetFixedPayFreq();
		if ( period == 12 )
		{
			itsBoostedVarTerm = string("1Y");
		}
		else
		{
			char str[2];
			sprintf(str, "%d", period);
			itsBoostedVarTerm = string(str) + string("M");
		}
		
		if (swaption->GetStrikes())
		{
			ARM_ReferenceValue* boostedStepUpFixRate = (ARM_ReferenceValue *) swaption->GetStrikes()->Clone();
			SetBoostedFixRate(*boostedStepUpFixRate);
			delete boostedStepUpFixRate;
		}
		else 
		{	
			double fixRate = fixedLeg->GetFixedRate();
			ARM_Vector* vecDates = new ARM_Vector(1, startDate.GetJulian());
			ARM_Vector* vecBoostedFixRate = new ARM_Vector(1, fixRate);
			ARM_ReferenceValue boostedFixRate = ARM_ReferenceValue(vecDates, vecBoostedFixRate);
			SetBoostedFixRate(boostedFixRate);
		}
/*		if (fixedLeg->GetFixedRates())
		{
			ARM_ReferenceValue* boostedStepUpFixRate = (ARM_ReferenceValue *) fixedLeg->GetFixedRates()->Clone();
			SetBoostedFixRate(*boostedStepUpFixRate);
			delete boostedStepUpFixRate;
		}
		else 
		{	
			double fixRate = fixedLeg->GetFixedRate();
			ARM_Vector* vecDates = new ARM_Vector(1, startDate.GetJulian());
			ARM_Vector* vecBoostedFixRate = new ARM_Vector(1, fixRate);
			ARM_ReferenceValue boostedFixRate = ARM_ReferenceValue(vecDates, vecBoostedFixRate);
			SetBoostedFixRate(boostedFixRate);
		}
*/	
		ARM_ReferenceValue* cpnSpread = new ARM_ReferenceValue(0.0);
		SetCpnSpread(*cpnSpread);
		delete cpnSpread;
	
		SetBoostedIndexType(boostedIndexType);
		
		//BoostedVarTerm
		string boostedVarTerm;
		int intBoostedVarTerm = cpnPayIndex->GetTerm();
		if (intBoostedVarTerm == K_ANNUAL)
			boostedVarTerm = "1Y";
		else if (intBoostedVarTerm == K_SEMIANNUAL)
			boostedVarTerm = "6M";
		else if (intBoostedVarTerm == K_QUARTERLY)
			boostedVarTerm = "3M";
		else if (intBoostedVarTerm == K_MONTHLY)
			boostedVarTerm = "1M";
		SetBoostedVarTerm(boostedVarTerm);
		
		//BoostedResetGap
		int boostedResetGap = cpnPayIndex->GetResetGap();
		SetBoostedResetGap(boostedResetGap);
		
		//BoostedResetTiming
		int boostedResetTiming = cpnPayIndex->GetResetTiming();
		SetBoostedResetTiming(boostedResetTiming);

		//BoostedDayCount
		int boostedDayCount = cpnPayIndex->GetDayCount(); 
		SetBoostedDayCount(boostedDayCount);
		
		//BoostedAdjRule
		int boostedAdjRule = cpnPayIndex->GetFwdRule();
		SetBoostedAdjRule(boostedAdjRule);
		
		//BoostedIntRule
		int boostedIntRule = cpnPayIndex->GetIntRule();
		SetBoostedIntRule(boostedIntRule);
		
		//CpnBarrierDown
		ARM_ReferenceValue* cpnBarrierDown = new ARM_ReferenceValue(-100.0);
		SetCpnBarrierDown(*cpnBarrierDown);
		delete cpnBarrierDown;
		
		//CpnBarrierUp
		ARM_ReferenceValue* cpnBarrierUp = new ARM_ReferenceValue(100.0);
		SetCpnBarrierUp(*cpnBarrierUp);
		delete cpnBarrierUp;
				
		//CpnResetFreq
		int cpnResetFreq = cpnRefIndex->GetResetFrequency();
		SetCpnResetFreq(cpnResetFreq);
		
		//CpnResetTiming
		int cpnResetTiming = cpnRefIndex->GetResetTiming();
		SetCpnResetTiming(cpnResetTiming);
		
		//NB: One of the two references indexes is not used.
		//RefIndexType
		int refIndexType;
		if (cpnRefIndex->IsLiborIndex())
			refIndexType = Libor;
		else
		{
			if (cpnRefIndex->IsCMSIndex())
				refIndexType = Cms;
		}
		SetRefIndexType1(refIndexType);
		
		//RefTerm
		string refTerm; 
		int intTerm = cpnRefIndex->GetTerm();
		if (intTerm == K_ANNUAL)
			refTerm = "1Y";
		else if (intTerm == K_SEMIANNUAL)
			refTerm = "6M";
		else if (intTerm == K_QUARTERLY)
			refTerm = "3M";
		else if (intTerm == K_MONTHLY)
			refTerm = "1M";
		SetRefTerm1(refTerm);
		
		//RefDayCount
		int refDayCount = cpnRefIndex->GetDayCount();
		SetRefDayCount1(refDayCount);
		
		SetRefCoeff1(1);
	
		//CallCal
		string callCal = cpnResetCal;
		SetCallCal(callCal);
		char* tempCallCal = (char *) callCal.c_str();
		
		//CallFees
		//Call notice : no need when we are in Summit case
		SetCallNotice(9999);

		//CallFreq
		int callFreq = cpnPayFreq; 
		SetCallFreq(callFreq);
	}
	catch (Exception& x)
	{
	    x.DebugPrint();

		throw x;
	}
	catch (...)
	{
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in CraCalculator::GenerateProductDescription(swaption)" );
	}
}

void ARM_CRACalculator::SetSwaption(ARM_Swaption* swaption)
{
	if (itsSwaption)
	{
		delete itsSwaption;
		itsSwaption = NULL;
	}
	if (swaption)
	{
		ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();
		if (swaption->GetStartDate() >= asOfDate)
		{
			itsSwaption = (ARM_Swaption*) swaption->Clone();
		}
		else
		{
			// keep only forward flows !!!
			ARM_SwapLeg* floatLeg = swaption->GetFloatLeg();
			ARM_SwapLeg* fixedLeg = swaption->GetFixedLeg();
			
			ARM_Date newStartDate;
			int i=0, j=0;
			// WARNING : if floatLeg frequency <> fixedLeg frequency, we may create a stub !
			if (floatLeg->GetIRIndex()->GetResetFrequency() >= fixedLeg->GetIRIndex()->GetPayFrequency())
			{
				while (asOfDate.GetJulian() > floatLeg->GetResetDates()->Elt(i))
				{
					i++;
					if (i >= floatLeg->GetResetDates()->size())
						ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : last floating leg period already fixed. Swaption should have at least one forward period.");
				}
				newStartDate = ARM_Date(floatLeg->GetFlowStartDates()->Elt(i));
			}
			else
			{
				while (asOfDate.GetJulian() > fixedLeg->GetResetDates()->Elt(j))
				{
					j++;
					if (j >= fixedLeg->GetResetDates()->size())
						ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : last fixed leg period already fixed. Swaption should have at least one forward period.");
				}
				newStartDate = ARM_Date(fixedLeg->GetFlowStartDates()->Elt(j));
			}

			if (swaption->GetExerciseStyle()->GetExerciseType() == K_EUROPEAN)
			{
				itsSwaption = new ARM_Swaption(newStartDate,
											   floatLeg->GetEndDateNA(),
											   swaption->GetOptionType(),
											   K_EUROPEAN,
											   swaption->GetStrike(),
											   swaption->GetMaturity(),
											   floatLeg->GetIRIndex()->GetIndexType(),
											   floatLeg->GetSpread(),
											   swaption->GetSwapYearTerm(),
											   floatLeg->GetIRIndex()->GetResetFrequency(),
											   floatLeg->GetIRIndex()->GetPayFrequency(),
											   swaption->GetCurrencyUnit());
			}
			else //BERMUDA
			{
				//Exercise Dates
				ARM_ExerciseStyle* newExerStyle = NULL;
				ARM_Vector* newExerDates = new ARM_Vector();
				ARM_Vector* newEndDates = new ARM_Vector();
				for (int i=0; i<swaption->GetExerciseStyle()->GetExerciseStartDates()->size(); i++)
				{
					if (asOfDate < ARM_Date(swaption->GetExerciseStyle()->GetExerciseStartDates()->Elt(i)))
					{
						newExerDates->push_back(swaption->GetExerciseStyle()->GetExerciseStartDates()->Elt(i));
						newEndDates->push_back(swaption->GetExerciseStyle()->GetExerciseEndDates()->Elt(i));
					}
				}
				newExerStyle = new ARM_ExerciseStyle(newExerDates,newEndDates);

				ARM_SwapLeg* fixLeg = new ARM_SwapLeg(newStartDate, 
													 fixedLeg->GetEndDateNA(), 
													 fixedLeg->GetFixedRate(), 
													 swaption->GetOptionType(),
													 fixedLeg->GetPaymentFreq(), 
													 fixedLeg->GetDayCount(), 
													 fixedLeg->GetDecompFreq(),
													 fixedLeg->GetIRIndex()->GetPayTiming(),
													 fixedLeg->GetIRIndex()->GetIntRule(),
													 fixedLeg->GetStubMeth(),
													 fixedLeg->GetCurrencyUnit(),
													 fixedLeg->GetPayCalName(),
													 fixedLeg->GetNxFlag(),
													 NULL, //ref date
													 1);   //adjust start

				ARM_SwapLeg* varLeg = new ARM_SwapLeg(newStartDate,
													 floatLeg->GetEndDateNA(),
													 floatLeg->GetIRIndex()->GetIndexType(),
													 swaption->GetOptionType(),
													 floatLeg->GetSpread(),
													 floatLeg->GetIRIndex()->GetResetFrequency(),
													 floatLeg->GetIRIndex()->GetPayFrequency(),
													 floatLeg->GetIRIndex()->GetResetTiming(),
													 floatLeg->GetIRIndex()->GetPayTiming(),
													 floatLeg->GetCurrencyUnit(),
													 floatLeg->GetIRIndex()->GetIntRule(),
													 floatLeg->GetIRIndex()->GetResetGap(),
													 floatLeg->GetResetCalName(),
													 floatLeg->GetPayCalName(),
													 floatLeg->GetDecompPricingFlag(),
													 floatLeg->GetNxFlag(),
													 floatLeg->GetStubMeth(),
													 NULL, //ref date
													 1,    //adjust start
													 floatLeg->GetDayCount(),
													 floatLeg->GetIRIndex()->GetFwdRule(),
													 floatLeg->GetIRIndex()->GetPayGap());

				if (floatLeg->GetVariableSpread())
				{
					varLeg->SetVariableSpread(floatLeg->GetVariableSpread());
				}

				ARM_Swap* swap = new ARM_Swap(fixLeg, varLeg);
			
				itsSwaption = new ARM_Swaption(swap, 
											   swaption->GetOptionType(),
											   newExerStyle,
											   swaption->GetStrikes(),
											   0);
				delete newExerStyle;
				delete fixLeg;
				delete varLeg;
				delete swap;
			}

			itsSwaption->SetAmount(floatLeg->GetAmount());
			itsSwaption->SetPorS(swaption->GetPorS());

#ifdef _DEBUG
			FILE* file = fopen("C:\\swaption.txt", "w");
			itsSwaption->View("", file);
			fclose(file);
#endif
		}
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRACalculator
///	Routine: CreateCalibMethod
///	Returns: void
///	Action : creates the calib method
/////////////////////////////////////////////////////////////////
void ARM_CRACalculator::CreateCalibMethods()
{
    try
	{
		ARM_ModelParams* Modelparams = &*GetPricingModel()->GetModelParams();
 
		ARM_ModelParamVector modelParamsVect;
		
		modelParamsVect = Modelparams->GetModelParams();

		///Create Calibmethods
		ARM_ModelParamVector calibParamSigmas(1, modelParamsVect[ARM_ModelParamType::Volatility]);
		ARM_CalibMethod* sigmaCalibMethod = new ARM_CalibMethod(itsSigmaPf,
													  calibParamSigmas,
													  ARM_CalibMethodType::Bootstrap1D, 
													  100, 
													  ARM_CalibrationTarget::PriceTarget,
													  NULL, 
													  NULL, 
													  false);

		if (itsToCalibrateBeta)
		{
			ARM_ModelParamVector calibParamBetas(1, modelParamsVect[ARM_ModelParamType::Beta]);        
    
			ARM_CalibMethod* calibMethodBeta = new ARM_CalibMethod(itsBetaPf,
									  calibParamBetas,
									  ARM_CalibMethodType::Optimize, 
									  100, 
									  ARM_CalibrationTarget::PriceTarget,
									  NULL, 
									  NULL,
									  false);
        
			if (itsToCalibrateMeanReversion)
			{
			   ARM_ModelParamVector calibParamMeanRevs(1, modelParamsVect[ARM_ModelParamType::MeanReversion]);
       
			   ARM_CalibMethod* calibMethodMR = new ARM_CalibMethod(itsMrsPf,
									   calibParamMeanRevs,
									   ARM_CalibMethodType::Optimize1D, 
									   100, 
									   ARM_CalibrationTarget::PriceTarget,
									   NULL, 
									   NULL,
									   false); 

			   SetCalibMethod(ARM_CalibMethodPtr(calibMethodMR));
			   GetCalibMethod()->SetlinkedMethod(calibMethodBeta);
			   GetCalibMethod()->GetlinkedMethod()->SetlinkedMethod(sigmaCalibMethod);
			}
			else
			{
				SetCalibMethod(ARM_CalibMethodPtr(calibMethodBeta));
				GetCalibMethod()->SetlinkedMethod(sigmaCalibMethod);
			}
		}
		else if (itsToCalibrateMeanReversion)
		{
		   ARM_ModelParamVector calibParamMeanRevs(1, modelParamsVect[1]);
       
		   ARM_CalibMethod* calibMethodMR = new ARM_CalibMethod(itsMrsPf,
								   calibParamMeanRevs,
								   ARM_CalibMethodType::Optimize1D, 
								   100, 
								   ARM_CalibrationTarget::PriceTarget,
								   NULL, 
								   NULL,
								   false); 
			
		   SetCalibMethod(ARM_CalibMethodPtr(calibMethodMR));
		   GetCalibMethod()->SetlinkedMethod(sigmaCalibMethod);
		}
		else
		{
			SetCalibMethod(ARM_CalibMethodPtr(sigmaCalibMethod));
		}
	}
	catch (Exception& x)
	{
	    x.DebugPrint();

		throw x;
	}
	catch (...)
	{
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in CraCalculator::CreateCalibMethods" );
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRACalculator
///	Routine: GetOSWCalibMethod
///	Returns: ARM_CalibMethod
///	Action : get calibration method for diagonal swaptions
/////////////////////////////////////////////////////////////////
ARM_CalibMethod* ARM_CRACalculator::GetOSWCalibMethod(void) const
{
	if( GetCalibMethod().IsNull() )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : No CalibMethod");

	return &(*GetCalibMethod());
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRACalculator
///	Routine: GetOSWPortfolio
///	Returns: ARM_Portfolio
///	Action : get the swaption calibration portfolio
/////////////////////////////////////////////////////////////////
const ARM_StdPortfolioPtr ARM_CRACalculator::GetOSWPortfolio(void) const
{
	if( GetCalibMethod()->GetPortfolio().IsNull() )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : No calibration porfolio");

	return GetCalibMethod()->GetPortfolio();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRACalculator
///	Routine: CreateAndSetCalibration
///	Returns: void
///	Action : create the calibration methods
/////////////////////////////////////////////////////////////////
void ARM_CRACalculator::CreateAndSetCalibration()
{
	try
	{	
	    CreateCalibMethods();
    }
	catch (Exception& x)
	{
	    x.DebugPrint();

		throw x;
	}
	catch (...)
	{
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in CraCalculator::CreateAndSetCalibration" );
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRACalculator
///	Routine: GetIndex
///	Returns: the index type
///	Action : get the ARM index type of the coupon index
/////////////////////////////////////////////////////////////////
ARM_IRIndex ARM_CRACalculator::GetIndex()
{
	if (skipCalib)
	{
		string liborTypeName(string(GetCurrencyUnit()->GetCcyName()) == string("EUR") ? "EURIBOR" : "LIBOR");

		string cpnIndexTerm = itsRefTerm1;

		if ( cpnIndexTerm == "12M" )
			cpnIndexTerm = "1Y"; // because we cant convert EURIBOR1Y and not EURIBOR12M

		int cpnIndexFreq=ARM_ArgConv_MatFrequency.GetNumber(cpnIndexTerm);
		int freq = CC_Max(cpnIndexFreq, itsCpnPayFreq);

		string term = ARM_ArgConvReverse_MatFrequency.GetString(freq);
		liborTypeName += term;

		ARM_IRIndex IRI(static_cast< ARM_INDEX_TYPE > (ARM_ArgConv_IndexType.GetNumber(liborTypeName)),
						K_DEF_FREQ,
						K_DEF_FREQ,
						GetCurrencyUnit());

		return(IRI);
	}
	else
	{
		ARM_CapFloor* capFloor = (ARM_CapFloor *) itsSigmaPf->GetAsset(0);

		return(*(capFloor->GetSwapLeg()->GetIRIndex()));
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRACalculator
///	Routine: CreateAndSetModel
///	Returns: void
///	Action : create the encapsulated model
/////////////////////////////////////////////////////////////////
void ARM_CRACalculator::CreateAndSetModel()
{
	try
	{
		ARM_ModelParamVector ModelParams;

		ARM_Model* capModel = dynamic_cast< ARM_Model* >( GetMktDataManager()->GetData(GetKeys()[CfModelKey]) );
		if (dynamic_cast<ARM_BSModel*>(capModel))
		  capModel = dynamic_cast<ARM_BSModel*>(capModel);
		else if (dynamic_cast<ARM_BSSmiledModel*>(capModel))
		  capModel = dynamic_cast<ARM_BSSmiledModel*>(capModel);

		ARM_CurveModelParam calibParamSigma;
		ARM_CurveModelParam calibParamMrs;
		ARM_CurveModelParam calibParamBeta;

		if (skipCalib) // for tests
		{
			std::vector<double> breakPointTimes(1, 0.0);
			std::vector<double> values(1, 0.0 ); // MeanRev = 0.0
			calibParamMrs = ARM_CurveModelParam(ARM_ModelParamType::MeanReversion, 
																		 &values, 
																		 &breakPointTimes);

			ModelParams.push_back(&calibParamMrs);

			values[0] = 1.0; // Beta = 100%
			calibParamBeta = ARM_CurveModelParam(ARM_ModelParamType::Beta, 
																		  &values, 
																		  &breakPointTimes); 

			ModelParams.push_back(&calibParamBeta);

			values[0] = 0.20; // Sigma = 20%
			calibParamSigma = ARM_CurveModelParam(ARM_ModelParamType::Volatility, 
																		   &values, 
																		   &breakPointTimes); 

			ModelParams.push_back(&calibParamSigma);
		}
		else
        {
            double asOfDate = capModel->GetStartDate().GetJulian();

            ARM_Vector* SigmaResetDates = !(itsSigmaPf.IsNull()) 
                                          ? itsSigmaPf->GetResetDates() 
                                          : NULL;
            
            ARM_Vector* BetaResetDates  = !(itsBetaPf.IsNull()) 
                                          ? itsBetaPf->GetResetDates() 
                                          : SigmaResetDates ? (ARM_Vector *) SigmaResetDates->Clone() : NULL;

            // use auto_ptr for exception safety!
            CC_NS(std,auto_ptr)<ARM_Vector> autoptrSigmaResetDates(SigmaResetDates);
            CC_NS(std,auto_ptr)<ARM_Vector> autoptrBetaResetDates(BetaResetDates);

            // Create model params
			size_t size_sig = itsSigmaPf->size();
			std::vector<double> sigmatimes(size_sig); 
			std::vector<double> sigmavalues(size_sig);

			size_t i;
    
			for (i = 0; i < size_sig; ++i)
				sigmatimes[i] = (*SigmaResetDates)[i] - asOfDate;

			for (i = 0; i < size_sig; ++i)
				sigmavalues[i] = itsSigmaCurve->Elt(i)/100.0;

			std::vector<double> lowerbound(size_sig, 1.0e-5);
			std::vector<double> upperbound(size_sig, UPPER_INFINITE_BOUND);

			calibParamSigma = ARM_CurveModelParam(ARM_ModelParamType::Volatility, 
												&sigmavalues, 
												&sigmatimes,
												"SIGMA", 
												"LINEAR",
												&lowerbound, 
												&upperbound);

			ModelParams.push_back(&calibParamSigma);

            // To create MRS params
			calibParamMrs = ARM_CurveModelParam(ARM_ModelParamType::MeanReversion, 
												  &std::vector<double>(1, (itsMeanRevMin + itsMeanRevMax)/2), 
												  &std::vector<double>(1, 0.0), 
												  "MEANREV", 
												  "STEPUPRIGHT",
												  &std::vector<double>(1, itsMeanRevMin), 
												  &std::vector<double>(1, itsMeanRevMax));

			ModelParams.push_back(&calibParamMrs);

            // To create Beta params
			size_t size_beta = BetaResetDates->GetSize();
			std::vector<double> betatimes(size_beta); 
			std::vector<double> betavalues(size_beta); 

			for (i = 0; i < size_beta; ++i)
				betatimes[i] = (*BetaResetDates)[i] - asOfDate;

			for (i = 0; i < size_beta; ++i)
				betavalues[i] = itsBetaOrShift->Elt(i);

			calibParamBeta = ARM_CurveModelParam(ARM_ModelParamType::Beta, 
											   &betavalues, 
											   &betatimes, 
											   "BETA", 
											   "LINEAR",
											   &std::vector<double>(size_beta, itsBetaMin),
											   &std::vector<double>(size_beta, itsBetaMax)); 

			ModelParams.push_back(&calibParamBeta);
        }

		// CreateModel

		ARM_IRIndex IRIndex = GetIndex();

		long nbFactors = 1; 

		int volType = itsVolType;

        // use the factory class to get the model params
		ARM_ModelParamsSFRM* pSFRMModelParams = 
			         ARM_ModelParamsSFRMFactory.Instance()->CreateModelParamsSFRM(ModelParams, 
                                                                                  &IRIndex, 
                                                                                  nbFactors, 
                                                                                  volType);

		/// use auto_ptr for exception safety!
		ARM_ZeroCurve* curve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));
	    CC_NS(std,auto_ptr)<ARM_ModelParamsSFRM> SFRMmodelParams(pSFRMModelParams);

        ARM_PricingModelPtr model;

        if (skipCalib)
        {
	       model = ARM_PricingModelPtr(static_cast< ARM_PricingModel* >(new ARM_SFRM(CreateClonedPtr(curve), 
                                       *pSFRMModelParams)));
        }
        else
        {
           ARM_SFRM* sfrmModel = new ARM_SFRM(CreateClonedPtr(curve), 
                                       *pSFRMModelParams,
                                       &*itsSigmaPf);
		   //sfrmModel->SetFixStartDate(itsStartDate);
		   model = ARM_PricingModelPtr(static_cast< ARM_PricingModel* >(sfrmModel));
        }

		// Tree Method
		int schedulerType = ARM_SchedulerBase::ConstantVariance;
		int samplerType = ARM_SamplerBase::NormalCentred;
		int truncatorType = ARM_TruncatorBase::StandardDeviation;
		int reconnectorType = ARM_ReconnectorBase::Mean;
		int smootherType = ARM_SmootherBase::DoNothing;

		std::vector<double> schedulerDatas(3);
		schedulerDatas[0] = itsNbSteps;
		schedulerDatas[1] = treeParams[1]; //should be = itsNbSteps/10

		std::vector<double> samplerDatas;
		std::vector<double> truncatorDatas(1,treeParams[2]);
		
		ARM_TreeBase* tree = ARM_TreeFactory.Instance()->CreateTreeND(SFRM_NB_FACTORS,
																	  schedulerType,
																	  schedulerDatas,
																	  samplerType,
																	  samplerDatas,
																	  truncatorType,
																	  truncatorDatas,
																	  false,
																	  reconnectorType,
																	  smootherType);

		model->SetNumMethod(ARM_NumMethodPtr(tree));

		// Numeraire: last neutral forward measure.
		ARM_NumerairePtr numeraire(ARM_NumeraireFactory.Instance()->CreateNumeraire(ARM_Numeraire::TerminalZc ) );
		model->SetNumeraire(numeraire);

		SetPricingModel(model);
	}

	catch (Exception& x)
	{
	    x.DebugPrint();

		throw x;
	}

	catch (...)
	{
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in CraCalculator::CreateAndSetModel" );
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRACalculator
///	Routine: SetReCalibFlags
///	Returns: void
///	Action : update the recalibration flags
///          update calibMethods (links)
/////////////////////////////////////////////////////////////////
void ARM_CRACalculator::SetReCalibFlags(int reCalibMrs, int reCalibBeta)
{
	if (itsReCalibMrs == 0 && reCalibMrs == 1)
		reCalibMrs = 0; // MRS has not been calibrated, so do not recalibrate !

	if (itsReCalibBeta == 0 && reCalibBeta == 1)
		reCalibBeta = 0; // Beta has not been calibrated, so do not recalibrate !

	// strategy has not changed => do nothing
	if ((itsReCalibMrs == reCalibMrs) && (itsReCalibBeta == reCalibBeta))
		return;

	ARM_CalibMethod* calibMethod = NULL;
	if ((itsReCalibMrs == reCalibMrs) && (itsReCalibBeta != reCalibBeta))
	{
		if (itsReCalibMrs == 0)
		{
			calibMethod = (ARM_CalibMethod*)(GetCalibMethod()->GetlinkedMethod()->Clone());
		}
		else
		{
			calibMethod = (ARM_CalibMethod*)(GetCalibMethod()->Clone());
			ARM_CalibMethod* linkedCalibMethod = (ARM_CalibMethod*)(GetCalibMethod()->GetlinkedMethod()->GetlinkedMethod()->Clone());
			calibMethod->SetlinkedMethod(linkedCalibMethod);
		}
	}
	if ((itsReCalibMrs != reCalibMrs) && (itsReCalibBeta == reCalibBeta))
	{
		calibMethod = (ARM_CalibMethod*)(GetCalibMethod()->GetlinkedMethod()->Clone());
	}
	if ((itsReCalibMrs != reCalibMrs) && (itsReCalibBeta != reCalibBeta))
	{
		calibMethod = (ARM_CalibMethod*)(GetCalibMethod()->GetlinkedMethod()->GetlinkedMethod()->Clone());
	}

	SetCalibMethod(ARM_CalibMethodPtr(static_cast< ARM_CalibMethod* >(calibMethod)));
	SetReCalibMrs(reCalibMrs);
	SetReCalibBeta(reCalibBeta);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRACalculator
///	Routine: UpdateModel
///	Returns: void
///	Action : update the model datas w.r.t. MarketDataManager
/////////////////////////////////////////////////////////////////
void ARM_CRACalculator::UpdateModel()
{
	ARM_ZeroCurve* cpnCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));

	ARM_SFRM* sfrmModel = dynamic_cast< ARM_SFRM* >(&*GetPricingModel());
	if( !sfrmModel)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : SFRM model is not a good type for updating");
  
	/// Update yield curves
	GetPricingModel()->SetZeroCurve( CreateClonedPtr( cpnCurve ) );
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRACalculator
///	Routine: UpdateCalibration
///	Returns: void
///	Action : update the calibration datas w.r.t. MarketDataManager
///          The context is an hedge ratios computation
/////////////////////////////////////////////////////////////////
void ARM_CRACalculator::UpdateCalibration(bool isUpdateStrike)
{
    try
	{
		/// Get the current market model for cap and swaption and its associated YC model
		ARM_Model* capModel = dynamic_cast< ARM_Model* >( GetMktDataManager()->GetData(GetKeys()[CfModelKey]) );
		if (dynamic_cast<ARM_BSModel*>(capModel))
		  capModel = dynamic_cast<ARM_BSModel*>(capModel);
		else if (dynamic_cast<ARM_BSSmiledModel*>(capModel))
		  capModel = dynamic_cast<ARM_BSSmiledModel*>(capModel);
		
		ARM_BSModel* swptModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswModelKey]) );

		/// Recompute market prices of caplets and swaptions		
		int i = 0;

		size_t size = itsSigmaPf->size();
		for (i = 0; i < size; ++i)
		{
			itsSigmaPf->GetAsset(i)->SetModelVariable(NULL);
			itsSigmaPf->GetAsset(i)->SetModel(capModel);

			double mktPrice = itsSigmaPf->GetAsset(i)->ComputePrice();
			double vega = itsSigmaPf->GetAsset(i)->ComputeSensitivity(K_VEGA);
			double curWeight = ( vega > 1.0e-5 ) ? 1.0 : 0.0;

			itsSigmaPf->SetWeight(curWeight, i);
			itsSigmaPf->SetPrecision(vega, i);
			itsSigmaPf->SetPrice(mktPrice, i);
		}
		
		//Beta portfolio
		if (itsReCalibBeta == 1)
		{
			size_t size = itsBetaPf->size();
			for (i = 0; i < size; ++i)
			{
				itsBetaPf->GetAsset(i)->SetModelVariable(NULL);
				itsBetaPf->GetAsset(i)->SetModel(capModel);

				double mktPrice = itsBetaPf->GetAsset(i)->ComputePrice();
				double vega = itsBetaPf->GetAsset(i)->ComputeSensitivity(K_VEGA);
				double curWeight = ( i == 0 ) ? 100 : 1.0;

				itsBetaPf->SetWeight(curWeight, i);
				itsBetaPf->SetPrecision(vega, i);
				itsBetaPf->SetPrice(mktPrice, i);
			}
		}
		//Mean reversion portfolio
		if (itsReCalibMrs == 1)
		{
			size_t size = itsMrsPf->size();
			for (i = 0; i < size; ++i)
			{
				itsMrsPf->GetAsset(i)->SetModelVariable(NULL);
				itsMrsPf->GetAsset(i)->SetModel(swptModel);

				double mktPrice = itsMrsPf->GetAsset(i)->ComputePrice();
				double vega = itsMrsPf->GetAsset(i)->ComputeSensitivity(K_VEGA);
				double curWeight = ( vega > 1.0e-5 ) ? 1.0 : 0.0;       

				itsMrsPf->SetWeight(curWeight, i);
				itsMrsPf->SetPrecision(vega, i);
				itsMrsPf->SetPrice(mktPrice, i);
			}
		}

		int szBeta = itsToCalibrateBeta ? itsBetaPf->GetSize() : itsSigmaPf->GetSize();
		double cstbeta = (itsBetaMax + itsBetaMin)/2.0; 
		itsBetaOrShift = new ARM_Vector(szBeta, cstbeta);	
	}
	catch (Exception& x)
	{
	    x.DebugPrint();

		throw x;
	}
	catch (...)
	{
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in CraCalculator::Update calibration" );
	}
}

ARM_PricingModel* ARM_CRACalculator::GetSFRMModel(void) const
{
	return &*(GetPricingModel());
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRACalculator
///	Routine: Calibrate
///	Returns: void
///	Action : calibrate the model 
/////////////////////////////////////////////////////////////////
void ARM_CRACalculator::Calibrate()
{
	try
	{
		if (!skipCalib)
		{
			itsReCalibMrs = itsReCalibMrs == -1 ? itsToCalibrateMeanReversion : itsReCalibMrs;

			if(itsReCalibMrs && (!itsToCalibrateMeanReversion))
				throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
									 "calibration flag are not valid : minMRS and MaxMrs should not be equal");
		
			itsReCalibBeta = itsReCalibBeta == -1 ? itsToCalibrateBeta : itsReCalibBeta;
			if(itsReCalibBeta && (!itsToCalibrateBeta))
				throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
									 "calibration flag are not valid : minBeta and MaxBeta should not be equal");

			///////////////////////////////////////////////////////////////////////////////////////////////
			///////////////////////Initialization of the Sigma Curve///////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////
			ARM_PricingModel* sfrmModel = GetSFRMModel();

			std::vector<double> mrs(1);

			mrs = (( (ARM_CurveModelParam&) (&*sfrmModel)->GetModelParams()->GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates());

			Pre_InitialiseFRMModel();

			ARM_FRMHWVol* FRMHWVol  = (ARM_FRMHWVol*)(itsFRMModel->GetFRMVol());

			ARM_ReferenceValue* sigmaRef  = FRMHWVol->GetInitSigmaCurve();
			size_t sizeref = sigmaRef->size();
			std::vector<double> sigmavalues(sizeref);
			int i;
			for(i = 0; i< sizeref; ++i)
				sigmavalues[i]  =(*(sigmaRef->GetDiscreteValues()))[i] / 100.;

			/// To set Volatility curve to Pricing Model and Calib Method
			( (ARM_CurveModelParam&) sfrmModel->GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility)).GetCurve()->SetOrdinates(sigmavalues);

			mrs = (( (ARM_CurveModelParam&) (&*sfrmModel)->GetModelParams()->GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates());

			///////////////////////////////////////////////////////////////////////////////////////////////
			///////////////////////Create CalibMethod if needed ///////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////		
			if (itsReCalibMrs == 1) 
			{
				std::vector<double> meanreversion(1, (itsMeanRevMin+itsMeanRevMax)/2);
				((ARM_CurveModelParam&) sfrmModel->GetModelParams()->GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->SetOrdinates(meanreversion);

				((ARM_CurveModelParam*) GetCalibMethod()->GetCalibParam(0))->GetCurve()->SetOrdinates(meanreversion);
				((ARM_CurveModelParam*) GetCalibMethod()->GetCalibParam(0))->GetInitialCurve()->SetOrdinates(meanreversion);
				GetCalibMethod()->SetPortfolio( ARM_StdPortfolioPtr((ARM_StdPortfolio*)(const_cast< ARM_StdPortfolio& >(*itsMrsPf).Clone())));

				if (itsReCalibBeta == 1) 
				{
					GetCalibMethod()->GetlinkedMethod()->GetlinkedMethod()->SetPortfolio( ARM_StdPortfolioPtr((ARM_StdPortfolio*)(const_cast< ARM_StdPortfolio& >(*itsSigmaPf).Clone())));
					((ARM_CurveModelParam*) &*GetCalibMethod()->GetlinkedMethod()->GetlinkedMethod()->GetCalibParam(0))->GetCurve()->SetOrdinates(sigmavalues);
					((ARM_CurveModelParam*) &*GetCalibMethod()->GetlinkedMethod()->GetlinkedMethod()->GetCalibParam(0))->GetInitialCurve()->SetOrdinates(sigmavalues);

					ARM_ReferenceValue* BetaCurve  = FRMHWVol->GetLogProbaCurve();
					size_t size_beta = BetaCurve->size();
				
					std::vector<double> betavalues(size_beta);
					for(i = 0; i< size_beta; ++i)
						betavalues[i]  =(*(BetaCurve->GetDiscreteValues()))[i];

					( (ARM_CurveModelParam&) sfrmModel->GetModelParams()->GetModelParam(ARM_ModelParamType::Beta)).GetCurve()->SetOrdinates(betavalues);
					((ARM_SFRM*)&*sfrmModel)->ConvertToShiftorBetaParam(*itsSigmaPf);

					((ARM_CurveModelParam*) GetCalibMethod()->GetlinkedMethod()->GetCalibParam(0))->GetCurve()->SetOrdinates(betavalues);
					((ARM_CurveModelParam*) GetCalibMethod()->GetlinkedMethod()->GetCalibParam(0))->GetInitialCurve()->SetOrdinates(betavalues);
					GetCalibMethod()->GetlinkedMethod()->SetPortfolio( ARM_StdPortfolioPtr((ARM_StdPortfolio*)(const_cast< ARM_StdPortfolio& >(*itsBetaPf).Clone())));
				}
				else
				{
					GetCalibMethod()->GetlinkedMethod()->SetPortfolio( ARM_StdPortfolioPtr((ARM_StdPortfolio*)(const_cast< ARM_StdPortfolio& >(*itsSigmaPf).Clone())));
					((ARM_CurveModelParam*) &*GetCalibMethod()->GetlinkedMethod()->GetCalibParam(0))->GetCurve()->SetOrdinates(sigmavalues);
					((ARM_CurveModelParam*) &*GetCalibMethod()->GetlinkedMethod()->GetCalibParam(0))->GetInitialCurve()->SetOrdinates(sigmavalues);
				}
			}

			else
			{
				if (itsReCalibBeta == 1) 
				{
					GetCalibMethod()->GetlinkedMethod()->SetPortfolio( ARM_StdPortfolioPtr((ARM_StdPortfolio*)(const_cast< ARM_StdPortfolio& >(*itsSigmaPf).Clone())));
					((ARM_CurveModelParam*) &*GetCalibMethod()->GetlinkedMethod()->GetCalibParam(0))->GetCurve()->SetOrdinates(sigmavalues);
					((ARM_CurveModelParam*) &*GetCalibMethod()->GetlinkedMethod()->GetCalibParam(0))->GetInitialCurve()->SetOrdinates(sigmavalues);
						
					ARM_ReferenceValue* BetaCurve  = FRMHWVol->GetLogProbaCurve();
					size_t size_beta = BetaCurve->size();
				
					std::vector<double> betavalues(size_beta);
					for(i = 0; i< size_beta; ++i)
						betavalues[i]  =(*(BetaCurve->GetDiscreteValues()))[i];

					( (ARM_CurveModelParam&) sfrmModel->GetModelParams()->GetModelParam(ARM_ModelParamType::Beta)).GetCurve()->SetOrdinates(betavalues);
					((ARM_SFRM*)&*sfrmModel)->ConvertToShiftorBetaParam(*itsSigmaPf);
				
					((ARM_CurveModelParam*) GetCalibMethod()->GetCalibParam(0))->GetCurve()->SetOrdinates(betavalues);
					((ARM_CurveModelParam*) GetCalibMethod()->GetCalibParam(0))->GetInitialCurve()->SetOrdinates(betavalues);
					GetCalibMethod()->SetPortfolio( ARM_StdPortfolioPtr((ARM_StdPortfolio*)(const_cast< ARM_StdPortfolio& >(*itsBetaPf).Clone())));
				}
				else
				{
					GetCalibMethod()->SetPortfolio( ARM_StdPortfolioPtr((ARM_StdPortfolio*)(const_cast< ARM_StdPortfolio& >(*itsSigmaPf).Clone())));
					((ARM_CurveModelParam*) &*GetCalibMethod()->GetCalibParam(0))->GetCurve()->SetOrdinates(sigmavalues);
					((ARM_CurveModelParam*) &*GetCalibMethod()->GetCalibParam(0))->GetInitialCurve()->SetOrdinates(sigmavalues);
				}
			}

			///////////////////////////////////////////////////////////////////////////////////////////////
			///////////////////////Calibrate Process///////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////	
			mrs = (( (ARM_CurveModelParam&) (&*sfrmModel)->GetModelParams()->GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates());

			GetCalibMethod()->Calibrate(&*sfrmModel);

			///////////////////////////////////////////////////////////////////////////////////////////////
			///////////////////////Update process//////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////		

			// set sigma calibrated curve
			sigmavalues = ( (ARM_CurveModelParam&) (&*sfrmModel)->GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility)).GetCurve()->GetOrdinates();
			for(i = 0; i< sizeref; ++i)
					(*(FRMHWVol->GetCurrentCurve()->GetDiscreteValues()))[i] = sigmavalues[i]* 100.;

			// set beta calibrated curve
			if (itsReCalibBeta == 1) 
			{
				std::vector<double> shift = ( (ARM_CurveModelParam&) (&*sfrmModel)->GetModelParams()->GetModelParam(ARM_ModelParamType::Shift)).GetCurve()->GetOrdinates();
				size_t size_shift=shift.size();
				if((*(FRMHWVol->GetMCurve()->GetDiscreteValues())).size() != size_shift)
					throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,"Inconsistency between shift curve in Kernel and GP");
				for(i = 0; i< size_shift; ++i)
					(*(FRMHWVol->GetMCurve()->GetDiscreteValues()))[i] = shift[i]* 100.;    

				std::vector<double> beta = ( (ARM_CurveModelParam&) (&*sfrmModel)->GetModelParams()->GetModelParam(ARM_ModelParamType::Beta)).GetCurve()->GetOrdinates();
				size_t size_beta=beta.size();
				if((*(FRMHWVol->GetLogProbaCurve()->GetDiscreteValues())).size() != size_beta)
					throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,"Inconsistency between shift curve in Kernel and GP");
				for(i = 0; i< size_beta; ++i)
					(*(FRMHWVol->GetLogProbaCurve()->GetDiscreteValues()))[i] = beta[i];

				SetReCalibBeta(1);
			}

			//set meanreversion
			if (itsReCalibMrs == 1) 
			{
				FRMHWVol->SetMeanRevParams((( (ARM_CurveModelParam&) (&*sfrmModel)->GetModelParams()->GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates())[0]);
				SetReCalibMrs(1);
			}

			for (int j = 0; j < sigmavalues.size(); ++j)    
					FRMHWVol->UpdateCurves(j);    

			itsFRMModel->SetFRMVol(FRMHWVol);
			itsFRMModel->UpDateErrors();

			bool isInitParam = true;
			if(isInitParam)
			{
				delete itsBetaOrShift;
				itsBetaOrShift   =(ARM_Vector*)itsFRMModel->GetFRMVol()->GetLogProbaCurve()->GetDiscreteValues()->Clone();
				itsMeanRev = itsFRMModel->GetFRMVol()->GetMeanRevParams();
			}
		}
	}
	catch (Exception& x)
	{
	    x.DebugPrint();

		throw x;
	}
	catch (...)
	{
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in CraCalculator::Calibrate" );
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRACalculator
///	Routine: Bump
///	Returns: 
///	Action : get & set the MRS param for SUMMIT interface
/////////////////////////////////////////////////////////////////
void ARM_CRACalculator::Bump(ARM_ZeroCurve* ZeroCurve, 
							ARM_VolCurve* CapVol,
							ARM_VolCurve* RhoCap,
							ARM_VolCurve* NuCap,
							ARM_VolCurve* SwoptVol,
							ARM_VolCurve* BetaCap,
							ARM_VolCurve* RhoSwopt,
							ARM_VolCurve* NuSwopt,
							ARM_VolCurve* BetaSwopt)
{
	double vSabrWeight = 0.5;
	int	vAlphaOrSigmaInput = 1;

	ARM_BSModel* CapModel = NULL;
	ARM_BSModel* SwoptModel = NULL;
	// Création du capModel
	if (RhoCap && NuCap)
	{
		int sabrCapFlag = K_SABR_ARITH;
		if (BetaCap)
		{
			if (BetaCap->IsEqualToOne() == false)
				sabrCapFlag = K_SABR_IMPLNVOL;
		}
	
		CapModel = new ARM_BSSmiledModel(ZeroCurve->GetAsOfDate(),
										0, //Spot
										ZeroCurve,
										ZeroCurve,
										CapVol,
										K_YIELD,
										RhoCap,
										NuCap,
										sabrCapFlag,
										BetaCap,
										vSabrWeight,
										vAlphaOrSigmaInput);
	}
	else
	{
		CapModel = new ARM_BSModel (ZeroCurve->GetAsOfDate(),
									0,
									ZeroCurve,
									ZeroCurve,
									CapVol,
									K_PRICE);
	}

	// Création du swoptModel
	if (RhoSwopt && NuSwopt)
	{
		int sabrSwoptFlag = K_SABR_ARITH;
		if (BetaSwopt)
		{
			if (BetaSwopt->IsEqualToOne() == false)
				sabrSwoptFlag = K_SABR_IMPLNVOL;
		}

		SwoptModel = new ARM_BSSmiledModel (ZeroCurve->GetAsOfDate(),
											0, //Spot
											ZeroCurve,
											ZeroCurve,
											SwoptVol,
											K_YIELD,
											RhoSwopt,
											NuSwopt,
											sabrSwoptFlag,
											BetaSwopt,
											vSabrWeight,
											vAlphaOrSigmaInput);
	}
	else
	{
		SwoptModel = new ARM_BSModel(ZeroCurve->GetAsOfDate(),
									 0,
									 ZeroCurve,
									 ZeroCurve,
									 SwoptVol,
									 K_PRICE);
	}
	
	ARM_MarketData_ManagerRep* mktDataManager = new ARM_MarketData_ManagerRep(ZeroCurve->GetAsOfDate());
	
	ARM_Currency* ccy = ZeroCurve->GetCurrencyUnit();
	string ccyName = ccy->GetCcyName();
	ARM_StringVector mdmKeys(4);
	mdmKeys[0] = string("YC_") + ccyName;
	mdmKeys[1] = string("OSWMOD_") + ccyName;
	mdmKeys[2] = string("CFMOD_") + ccyName;
	mdmKeys[3] = string("MRS_") + ccyName;
	
	mktDataManager->RegisterData(mdmKeys[0], ZeroCurve);
	mktDataManager->RegisterData(mdmKeys[1], SwoptModel);
	mktDataManager->RegisterData(mdmKeys[2], CapModel);

	std::vector<double> mrs(1, itsMeanRev );
	std::vector<double> breakPointTimes(1, 0.0);
	ARM_CurveModelParam* MeanRev = new ARM_CurveModelParam (ARM_ModelParamType::MeanReversion, 
															&mrs, 
															&breakPointTimes);
	
	mktDataManager->RegisterData(mdmKeys[3], MeanRev);

	// Market data have changed, so reset them to calculator...
	SetMktDataManager(ARM_MarketData_ManagerRepPtr(mktDataManager));
	// ... and update model and calibration
	Update();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRACalculator
///	Routine: GetMRS & SetMRS
///	Returns: 
///	Action : get & set the MRS param for SUMMIT interface
/////////////////////////////////////////////////////////////////
const ARM_ModelParam& ARM_CRACalculator::GetMRS() const
{
    return *(static_cast< ARM_ModelParam* >(GetMktDataManager()->GetData(GetKeys()[MrsKey])));
}

void ARM_CRACalculator::SetMRS(ARM_ModelParam* mrsParam)
{
    if(!mrsParam || mrsParam->GetType() != ARM_ModelParamType::MeanReversion)
		    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : an MRS Param is expected !" );

    GetMktDataManager()->RegisterData(GetKeys()[MrsKey],static_cast< ARM_Object* >(mrsParam));
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRACalculator
///	Routine: GetBeta & SetBeta
///	Returns: 
///	Action : get & set the Beta param for SUMMIT interface
/////////////////////////////////////////////////////////////////
const ARM_ModelParam& ARM_CRACalculator::GetBeta() const
{
    return *(static_cast< ARM_ModelParam* >(GetMktDataManager()->GetData(GetKeys()[BetaKey])));
}

void ARM_CRACalculator::SetBeta(ARM_ModelParam* betaParam)
{
    if(!betaParam || betaParam->GetType() != ARM_ModelParamType::Beta)
		    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : an Beta Param is expected !" );

    GetMktDataManager()->RegisterData(GetKeys()[BetaKey],static_cast< ARM_Object* >(betaParam));
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRACalculator
///	Routine: Price
///	Returns: a double
///	Action : price the CRA.
/////////////////////////////////////////////////////////////////
double ARM_CRACalculator::Price()
{
	try
	{
		CalibrateAndTimeIt();

		/// Price the implicit product according to internal flag
		ARM_GenSecurityPtr genSec = GetGenSecurity();

		size_t nbPricedColumns = genSec->GetDealDescription().GetPricedColumnNames().size();

		double price = 0.0;
		ARM_GenPricer* genPricer = NULL;
		genPricer = new ARM_GenPricer( &*genSec,&*GetPricingModel() );
		ARM_AutoCleaner<ARM_GenPricer> HoldGP(genPricer );
		genPricer->Price();

		// effective nb of columns to price :
		int colSize = MIN(NbProductsToPrice, itsProductsToPrice.size());

		for (size_t i(0); i < colSize; ++i)
		{
			if (itsProductsToPrice[i])
			{
				price	= genPricer->GetPricerInfo()->GetContents( CRAColNamesTable[ CRAProductToPriceColumns[i] ] ).GetData("Price").GetDouble();
						
				if (i == BermudaPrice)
				{
					itsBermuda1Price = price;
				}
				else if (i == Bermuda2Price)
				{
					itsBermuda2Price = price;
				}
				else if (i == FundingPrice)
				{
					itsFundingPrice = price;
				}
				else if (i == CorridorPrice)
				{
					itsCorridorPrice = price;
				}
				else if (i == ExoticSwapPrice)
				{
					itsExoticSwapPrice = price;
				}
				else if (i == FwdCorridorPrice)
				{
					itsFwdCorridorPrice = price;
				}
				else if (i == FwdFundingPrice)
				{
					itsFwdFundingPrice = price;
				}
			}
			itsHasBeenPriced = true;
		}

		price = itsBermuda1Price;
		return price;
	}
	catch (Exception& x)
	{
	    x.DebugPrint();

		throw x;
	}
	catch (...)
	{
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in CraCalculator::Price" );
	}
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRACalculator
///	Routine: ComputePricingData
///	Returns: a ARM_MultiTypeDict
///	Action : get pricing data of the calculaltor
/////////////////////////////////////////////////////////////////
void ARM_CRACalculator::ComputePricingData() const
{
	if (!itsHasBeenPriced)
		const_cast<ARM_CRACalculator*>(this)->PriceAndTimeIt();

	GetPricingData()[ "Bermuda" ]			= itsBermuda1Price;
	GetPricingData()[ "Bermuda2" ]			= itsBermuda2Price;
	GetPricingData()[ "Funding" ]			= itsFundingPrice;		
	GetPricingData()[ "Corridor" ]			= itsCorridorPrice;		
	GetPricingData()[ "ExoticSwap" ]		= itsExoticSwapPrice;	
	GetPricingData()[ "FwdCorridor" ]		= itsFwdCorridorPrice;		
	GetPricingData()[ "FwdFunding" ]		= itsFwdFundingPrice;	
}

////////////////////////////////////////////////////
///	Class   : ARM_CRACalculator
///	Routines: toString
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
string ARM_CRACalculator::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream craData;

    return craData.str();
}

////////////////////////////////////////////////////
///	Class   : ARM_CRACalculator
///	Routines: GeneralDataToString
///	Returns :
///	Action  : Construct information string for
///           general product data
////////////////////////////////////////////////////
string ARM_CRACalculator::GeneralDataToString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream generalCraData;

	//General Datas
	generalCraData << "GENERAL DATAS"	<< endl;
	generalCraData << "Currency\t\t: "	<<  itsCcy.GetCcyName() << endl;
	generalCraData << "StartDate\t\t: " <<  itsStartDate.toString() << endl;
    generalCraData << "EndDate\t\t: "	<< itsEndDate.toString() << endl;
    generalCraData << "Pay/Rec\t\t: "	<< ARM_ParamView::GetMappingName(S_RECEIVE_PAY, itsPayRec) << endl;
	generalCraData << "\n";
	
	//Call Datas
	generalCraData << "CALL DATAS" << endl;
	generalCraData << "\nCall Frequency\t\t: "	<< ARM_ArgConvReverse_StdFrequency.GetString(itsCallFreq) << endl;
   	generalCraData << "Call Notice\t\t:  "		<< itsCallNotice << "\n";
	generalCraData << "Call Calendar\t\t: "		<< itsCallCal << "\n";
	generalCraData << "\n";

	//Fund Datas
	generalCraData << "FUND DATAS" << endl;
    generalCraData << "Fund Frequency\t\t: "	<< ARM_ArgConvReverse_StdFrequency.GetString(itsFundFreq) << endl;
	generalCraData << "Fund Day Count\t\t: "	<< ARM_ArgConvReverse_DayCount.GetString(itsFundDayCount)  << endl;
	generalCraData << "\n";

	//Coupon Datas
	generalCraData << "EXOTIC COUPON DATAS"	<< endl;
    generalCraData << "Cpn Pay Freq\t\t: "			<< ARM_ArgConvReverse_StdFrequency.GetString(itsCpnPayFreq) << endl;
   	generalCraData << "Cpn Reset Calendar\t\t: "	<< itsCpnResetCal << "\n";
	generalCraData << "Cpn Pay Calendar\t\t: "		<< itsCpnPayCal << "\n";

	generalCraData << "Boosted Index Type\t\t: "	<< ARM_ArgConvReverse_IndexClass.GetString(itsBoostedIndexType) << endl;
	
	if (itsBoostedIndexType != K_FIXED)
	{
		generalCraData << "Boosted Var Term\t\t: "	<< itsBoostedVarTerm << endl;
		generalCraData << "Boosted Reset Gap\t\t: " << itsBoostedResetGap << endl;
		generalCraData << "Boosted Reset Timing\t\t: " << ARM_ArgConvReverse_Timing.GetString(itsBoostedResetTiming) << endl;
		generalCraData << "Boosted Day Count\t\t: " << ARM_ArgConvReverse_DayCount.GetString(itsBoostedDayCount)  << endl;
		generalCraData << "Boosted Adj Rule\t\t: "	<< ARM_ArgConvReverse_FwdRules.GetString(itsBoostedAdjRule)  << endl;
		generalCraData << "Boosted Rule\t\t: "		<< ARM_ArgConvReverse_InterestRules.GetString(itsBoostedIntRule)  << endl;
	}

	generalCraData << "Cpn Reset Freq\t\t: "		<< ARM_ArgConvReverse_StdFrequency.GetString(itsCpnResetFreq) << endl;
    generalCraData << "Cpn Reset Timing\t\t: "		<< ARM_ArgConvReverse_Timing.GetString(itsCpnResetTiming) << endl;
	generalCraData << "Cpn Reset Gap\t\t: "			<< itsCpnResetGap << "\n";

	generalCraData << "Ref Index Type 1\t\t: "		<< ARM_ArgConvReverse_IndexClass.GetString(itsRefIndexType1) << endl;
	generalCraData << "Ref Term 1\t\t: "			<< itsRefTerm1 << "\n";
	generalCraData << "Ref Day Count 1\t\t: "		<< ARM_ArgConvReverse_DayCount.GetString(itsRefDayCount1)  << endl;
	generalCraData << "Ref Coeff\t\t: "				<< itsRefCoeff1 << "\n";

	generalCraData << endl;

    return generalCraData.str();
}

////////////////////////////////////////////////////
///	Class   : ARM_CRACalculator
///	Routines: DealDesDataToString
///	Returns :
///	Action  : Construct information string for
/// deal description data
////////////////////////////////////////////////////
string ARM_CRACalculator::DealDesDataToString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream dealDesCraData;

	double asOf = GetMktDataManager()->GetAsOfDate().GetJulian();
	//Notional Profile
	dealDesCraData << " ===========Notional Profile=============" << endl;
	ARM_GP_CurvePtr notional(RefValueToCurve(itsNotional, asOf));
	dealDesCraData << notional->toString();

	//Call Fees Profile
	dealDesCraData << " ===========Call Fees Profile=============" << endl;
	ARM_GP_CurvePtr callFees(RefValueToCurve(itsCallFees, asOf));
	dealDesCraData << callFees->toString();

	//Fund Spread Profile
	dealDesCraData << " ===========Fund Spread Profile=============" << endl;
	ARM_GP_CurvePtr fundSpread(RefValueToCurve(itsFundSpread, asOf));
	dealDesCraData << fundSpread->toString();

	//Cpn Spread Profile
	// removed since cpnSpread seems to be never used
	if (false)///if (itsBoostedIndexType != K_FIXED)
	{
		dealDesCraData << " ===========Cpn Spread Profile=============" << endl;
		ARM_GP_CurvePtr cpnSpread(RefValueToCurve(itsCpnSpread, asOf));
		dealDesCraData << cpnSpread->toString();
	}

	//Boosted Fix Rate Profile
	dealDesCraData << " ===========Boosted Fix Rate Profile=============" << endl;
	ARM_GP_CurvePtr boostedFixRate(RefValueToCurve(itsBoostedFixRate, asOf));
	dealDesCraData << boostedFixRate->toString();

	//Cpn Barrier Down
	dealDesCraData << " ===========Cpn Barrier Down Profile=============" << endl;
	ARM_GP_CurvePtr cpnBarrierDown(RefValueToCurve(itsCpnBarrierDown, asOf));
	dealDesCraData << cpnBarrierDown->toString();

	//Cpn Barrier Up
	dealDesCraData << " ===========Cpn Barrier Up Profile=============" << endl;
	ARM_GP_CurvePtr cpnBarrierUp(RefValueToCurve(itsCpnBarrierUp, asOf));
	dealDesCraData << cpnBarrierUp->toString();

    return dealDesCraData.str();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRACalculator
///	Routine: View
///	Returns: 
///	Action : .
/////////////////////////////////////////////////////////////////
void ARM_CRACalculator::View(char* id, FILE* ficOut) const
{
	FILE* fOut;
    char fOutName[200];
	
    if ( ficOut == NULL )
    {
		ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
		(void) unlink(fOutName);
		fOut = fopen(fOutName, "w");
    }
    else
		fOut = ficOut;

    /// Callable Range Accrual Calculator specific datas viewing
    fprintf(fOut,"\n\n ===================================================== \n");
	fprintf(fOut,"\n\n ========> CALLABLE RANGE ACCRUAL CALCULATOR <======== \n");
    fprintf(fOut,"\n\n ===================================================== \n");

	fprintf(fOut,"%s",GeneralDataToString().c_str());
	fprintf(fOut,"%s",DealDesDataToString().c_str());

	//Calib Method / Portfolio
	CC_Ostringstream calibPf;
	fprintf(fOut,"\n\n ============>    CALIB - PORTFOLIO      <=========== \n");
	calibPf << GetCalibMethod()->GetPortfolio()->toString() << "\n";
	calibPf << GetCalibMethod()->toString() << "\n";

	if (GetCalibMethod()->GetlinkedMethod())
	{
		calibPf << GetCalibMethod()->GetlinkedMethod()->GetPortfolio()->toString() << "\n";
		calibPf << GetCalibMethod()->GetlinkedMethod()->toString() << "\n";
	    if ( GetCalibMethod()->GetlinkedMethod()->GetlinkedMethod())
	    {
		    calibPf << GetCalibMethod()->GetlinkedMethod()->GetlinkedMethod()->GetPortfolio()->toString() << "\n";
		    calibPf << GetCalibMethod()->GetlinkedMethod()->GetlinkedMethod()->toString() << "\n";
	    }
	}
	fprintf(fOut,"%s",calibPf.str().c_str());
	fprintf(fOut,"\n\n ============>  EOF CALIB - PORTFOLIO <=========== \n");
	
	fprintf(fOut,"\n\nRecalibrate MRS flag:%d\n", itsReCalibMrs);
	fprintf(fOut,"Recalibrate Beta flag:%d\n", itsReCalibBeta);

	CC_Ostringstream genData;
	fprintf(fOut,"\n\n ============>    GEN CALCULATOR    <=========== \n");
	genData << "\n\nCommon Calculator Part\n" << ARM_GenCalculator::toString() << "\n\n";
	fprintf(fOut,"%s",genData.str().c_str());
	fprintf(fOut,"\n\n ============>  EOF GEN CALCULATOR  <=========== \n");

	fprintf(fOut,"\n\n ============>    EQUIVALENT OPTION PORTFOLIO    <=========== \n");
	if (itsOptionPortfolio)
		itsOptionPortfolio->View(id, fOut);
	fprintf(fOut,"\n\n ============>  EOF EQUIVALENT OPTION PORTFOLIO  <=========== \n");

	fprintf(fOut,"\n\n ======END OF VIEW=============================== \n");

	if ( ficOut == NULL )
       fclose(fOut);
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
