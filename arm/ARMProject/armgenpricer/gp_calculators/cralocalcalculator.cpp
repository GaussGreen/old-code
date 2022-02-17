
/*!
 /** Copyright (c) CDC IXIS CM July 2003 Paris
 *	\file CRALocalCalculator.cpp
 *  \brief file for the Local CRA Calculator
 *	\author  H. BAKHTRI & M. ABDELMOUMNI & P. LAM
 *	\version 1.0
 *	\date August 2005
 */

/// this header comes first as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpcalculators/cralocalcalculator.h"
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
#include "gpinfra/surfacelistmodelparam.h"
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
#include "gpmodels/local_SLN_Model.h"
#include "gpmodels/local_SLN_ModelParams.h"
#include "gpmodels/MultiAssetsFactory.h"
#include "gpmodels/multiassets.h"

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

/// Call schedule 
const unsigned int CALL_CRA_SCHED	= 0;
const double NON_CALL_FEE			= 1.0e15;

const int SFRM_NB_FACTORS			= 1;

const double treeParams[]			= {3.0, 30.0, 7.0, 2.0, 0.0};

const string ARM_CRALocalCalculator::LocalFloorCRAColNamesTable [] =
{  
	"ResetDate",
	"StartDate",
	"MaturityDate",
	"Fees",
	"FundingLeg",
	"CorridorDown",
	"CorridorUp",
	"CorridorLeg",
	"ExoticSwap",
	"Option", 	
	"Bermuda",
	"Corridor",
	"Funding",
	"FwdCorridorUp",
	"FwdCorridorDown",
	"FwdCorridor",
	"FwdFunding",
};

const size_t ARM_CRALocalCalculator::LocalCapFloorCRAExerciseProbaNb=20;
const string ARM_CRALocalCalculator::LocalCapFloorCRALocalExerciseProbaNamesTable [] =
{
	"LocalProba1",
	"LocalProba2",
	"LocalProba3",
	"LocalProba4",
	"LocalProba5",
	"LocalProba6",
	"LocalProba7",
	"LocalProba8",
	"LocalProba9",
	"LocalProba10",
	"LocalProba11",
	"LocalProba12",
	"LocalProba13",
	"LocalProba14",
	"LocalProba15",
	"LocalProba16",
	"LocalProba17",
	"LocalProba18",
	"LocalProba19",
	"LocalProba20"
};
const string ARM_CRALocalCalculator::LocalCapFloorCRAExerciseProbaNamesTable [] =
{
	"Proba1",
	"Proba2",
	"Proba3",
	"Proba4",
	"Proba5",
	"Proba6",
	"Proba7",
	"Proba8",
	"Proba9",
	"Proba10",
	"Proba11",
	"Proba12",
	"Proba13",
	"Proba14",
	"Proba15",
	"Proba16",
	"Proba17",
	"Proba18",
	"Proba19",
	"Proba20"
};

const string ARM_CRALocalCalculator::LocalCapFloorCRAColNamesTable [] =
{  
	"ResetDate",
	"StartDate",
	"MaturityDate",
	"Fees",
	"FundingLeg",
	"CorridorLeg1",
	"CorridorLeg2",
	"CorridorLeg3",
	"CorridorLeg",
	"ExoticSwap",
	"Option", 
	"Bermuda",
	"Corridor",
	"Funding",
	"FwdCorridor",
	"FwdFunding",
	"ExerSwapRate",
	"Frontier",
	"NextResetDate",
	"IndicExer",

	LocalCapFloorCRALocalExerciseProbaNamesTable[0],
	LocalCapFloorCRAExerciseProbaNamesTable[0],

	LocalCapFloorCRALocalExerciseProbaNamesTable[1],
	LocalCapFloorCRAExerciseProbaNamesTable[1],

	LocalCapFloorCRALocalExerciseProbaNamesTable[2],
	LocalCapFloorCRAExerciseProbaNamesTable[2],

	LocalCapFloorCRALocalExerciseProbaNamesTable[3],
	LocalCapFloorCRAExerciseProbaNamesTable[3],

	LocalCapFloorCRALocalExerciseProbaNamesTable[4],
	LocalCapFloorCRAExerciseProbaNamesTable[4],

	LocalCapFloorCRALocalExerciseProbaNamesTable[5],
	LocalCapFloorCRAExerciseProbaNamesTable[5],

	LocalCapFloorCRALocalExerciseProbaNamesTable[6],
	LocalCapFloorCRAExerciseProbaNamesTable[6],

	LocalCapFloorCRALocalExerciseProbaNamesTable[7],
	LocalCapFloorCRAExerciseProbaNamesTable[7],

	LocalCapFloorCRALocalExerciseProbaNamesTable[8],
	LocalCapFloorCRAExerciseProbaNamesTable[8],

	LocalCapFloorCRALocalExerciseProbaNamesTable[9],
	LocalCapFloorCRAExerciseProbaNamesTable[9],

	LocalCapFloorCRALocalExerciseProbaNamesTable[10],
	LocalCapFloorCRAExerciseProbaNamesTable[10],

	LocalCapFloorCRALocalExerciseProbaNamesTable[11],
	LocalCapFloorCRAExerciseProbaNamesTable[11],

	LocalCapFloorCRALocalExerciseProbaNamesTable[12],
	LocalCapFloorCRAExerciseProbaNamesTable[12],

	LocalCapFloorCRALocalExerciseProbaNamesTable[13],
	LocalCapFloorCRAExerciseProbaNamesTable[13],

	LocalCapFloorCRALocalExerciseProbaNamesTable[14],
	LocalCapFloorCRAExerciseProbaNamesTable[14],

	LocalCapFloorCRALocalExerciseProbaNamesTable[15],
	LocalCapFloorCRAExerciseProbaNamesTable[15],

	LocalCapFloorCRALocalExerciseProbaNamesTable[16],
	LocalCapFloorCRAExerciseProbaNamesTable[16],

	LocalCapFloorCRALocalExerciseProbaNamesTable[17],
	LocalCapFloorCRAExerciseProbaNamesTable[17],

	LocalCapFloorCRALocalExerciseProbaNamesTable[18],
	LocalCapFloorCRAExerciseProbaNamesTable[18],

	LocalCapFloorCRALocalExerciseProbaNamesTable[19],
	LocalCapFloorCRAExerciseProbaNamesTable[19]
};

const int ARM_CRALocalCalculator::LocalFloorCRAProductToPriceColumns [] =
{
    localFBermuda,
	localFCorridor,
	localFFunding,
	localFFwdCorridor,
	localFFwdFunding,
};

const int ARM_CRALocalCalculator::LocalCapFloorCRAProductToPriceColumns [] =
{
	localCfBermuda,
	localCfCorridor,
	localCfFunding,
	localCfFwdCorridor,
	localCfFwdFunding,
	localCfCorridorLeg
};

////////////////////////////////////////////////////
///	Class  : ARM_CRALocalCalculator
///	Routine: CopyNoCleanUp
///	Returns: void
///	Action : Arguments copy
////////////////////////////////////////////////////
void ARM_CRALocalCalculator::CopyNoCleanUp(const ARM_CRALocalCalculator& rhs)
{
	ARM_CRACalculator::CopyNoCleanUp(rhs);

	itsLocalModelType		= rhs.itsLocalModelType;
	itsLocalResetFreq		= rhs.itsLocalResetFreq; 
	itsLocalPortfolio		= rhs.itsLocalPortfolio;
	itsIsExerciseProbas		= rhs.itsIsExerciseProbas;
	itsExerciseProbaOffset	= rhs.itsExerciseProbaOffset;
	itsExerciseProbas		= rhs.itsExerciseProbas;
}

////////////////////////////////////////////////////
///	Class  : ARM_CRALocalCalculator
///	Routine: CleanUp
///	Returns: void
///	Action : delete pointors
////////////////////////////////////////////////////
void ARM_CRALocalCalculator::CleanUp()
{
	ARM_CRACalculator::CleanUp();
}

ARM_CRALocalCalculator::ARM_CRALocalCalculator(const ARM_MarketData_ManagerRep& mktDataManager)
					   :ARM_CRACalculator(mktDataManager),
itsIsExerciseProbas(false),
itsExerciseProbaOffset(0),
itsExerciseProbas(std::vector<double>(1,0.0))
{
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRALocalCalculator /////////////////////////////
///	Routine: Constructor ////////////////////////////////////////
///	Returns: void ///////////////////////////////////////////////
///	Action : builds the object //////////////////////////////////
/////////////////////////////////////////////////////////////////
ARM_CRALocalCalculator::ARM_CRALocalCalculator( ARM_Currency& ccy,
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
												int localModelType,
												double meanRevMin,
												double meanRevMax,
												double betaMin,
												double betaMax,
												ARM_Vector* calibSecPFParams,
												int nbSteps,
												int flagToGenerateOSWATM,
												int localResetFreq,
												ARM_StringVector& mdmKeys,
												const ARM_MarketData_ManagerRep& mktDataManager,
												const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice,
												int reCalibMrs,
												int reCalibBeta,
												bool isStdCalib,
												bool isExerciseProbas,
												size_t exerciseProbaOffset,
												ARM_Currency* fundccy,
												const ARM_ReferenceValue& fundNominal )		  
: ARM_CRACalculator(mktDataManager),
itsLocalModelType(localModelType),
itsLocalResetFreq(localResetFreq),
itsLocalPortfolio(ARM_StdPortfolioPtr(NULL)),
itsIsExerciseProbas(isExerciseProbas),
itsExerciseProbaOffset(exerciseProbaOffset),
itsExerciseProbas(std::vector<double>(1,0.0))
{
	Set( ccy, startDate, endDate, payRec, notional, callFreq, callNotice,
	     callCal, callFees, fundFreq, fundSpread, fundDayCount, cpnSpread,
	     cpnPayFreq, cpnResetCal, cpnPayCal, boostedIndexType,
	     boostedFixRate, boostedVarTerm, boostedResetGap, boostedResetTiming,
	     boostedDayCount, boostedAdjRule, boostedIntRule,
	     cpnBarrierDown, cpnBarrierUp, cpnResetFreq, cpnResetTiming, cpnResetGap,
	     refIndexType, refTerm, refDayCount, refCoeff, 
	     meanRevMin, meanRevMax, betaMin, betaMax, nbSteps, calibSecPFParams,
	     flagToGenerateOSWATM, mdmKeys, productsToPrice, reCalibMrs, reCalibBeta, K_DIAG,isStdCalib);

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

	SetKeys(mdmKeys);

	CheckData();

	ARM_StringVector pricedColumns = PricedColumnNames();

	ARM_CstManagerPtr cstManagerPtr = CreateCstManager();

	Initialize();

	if (itsLocalModelType == LocalDown)
		CreateAndSetDealDescriptionAndTimeIt("LOCALFLOOR", pricedColumns, cstManagerPtr);
	else if (itsLocalModelType == LocalUp)
		CreateAndSetDealDescriptionAndTimeIt("LOCALCAP", pricedColumns, cstManagerPtr);
	else 
		CreateAndSetDealDescriptionAndTimeIt("LOCALCAPFLOOR", pricedColumns, cstManagerPtr);

	//SFRM Model + Local Model
	CreateAndSetModel();

	//SFRM Calib Methods
	CreateAndSetCalibration();
	
	//LOCAL CORRIDOR PORTFOLIO
	CreateCorridorLocalPortfolio();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRALocalCalculator /////////////////////
///	Routine: Summit Constructor from Option Portfolio ///////////
///	Returns: void ///////////////////////////////////////////////
///	Action : builds the object //////////////////////////////////
/////////////////////////////////////////////////////////////////
ARM_CRALocalCalculator::ARM_CRALocalCalculator( ARM_OptionPortfolio* optPortfolio,
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
												int localModelType,
												int localResetFreq,
												bool isStdCalib)
:ARM_CRACalculator(mktDataManager),
itsLocalResetFreq(localResetFreq),
itsLocalPortfolio(ARM_StdPortfolioPtr(NULL)),
itsIsExerciseProbas(false),
itsExerciseProbaOffset(0),
itsExerciseProbas(std::vector<double>(1,0.0))
{
	Set(optPortfolio, meanRevMin, meanRevMax, betaMin, betaMax,
		calibSecPFParams, nbSteps, flagToGenerateOSWATM,
		mdmKeys, productsToPrice, reCalibMrs, reCalibBeta, K_DIAG, isStdCalib);

	SetLocalModelType(localModelType);

	CheckCRAInputs();

	CheckData();

	GenerateProductDescription((ARM_OptionPortfolio *) optPortfolio);

	ARM_StringVector pricedColumns = PricedColumnNames();

	ARM_CstManagerPtr cstManagerPtr = CreateCstManager();

	Initialize();

	if (itsLocalModelType == LocalDown)
		CreateAndSetDealDescriptionAndTimeIt("LOCALFLOOR", pricedColumns, cstManagerPtr);
	else if (itsLocalModelType == LocalUp)
		CreateAndSetDealDescriptionAndTimeIt("LOCALCAP", pricedColumns, cstManagerPtr);
	else
		CreateAndSetDealDescriptionAndTimeIt("LOCALCAPFLOOR", pricedColumns, cstManagerPtr);

	CreateAndSetModel();

	CreateCalibMethods(); 
	
	CreateCorridorLocalPortfolio();
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRALocalCalculator 
///	Routine: Summit Constructor from Option Portfolio 
///	Returns: void
///	Action : builds the object
/// WARNING: this constructor is called ONLY from CRA Spread
///          to degenerate the CCSO into CRA
///          then no need to initialise the SFRM and local models
/////////////////////////////////////////////////////////////////
ARM_CRALocalCalculator::ARM_CRALocalCalculator( ARM_OptionPortfolio* optPortfolio,
												ARM_StringVector& mdmKeys,
												const ARM_MarketData_ManagerRep& mktDataManager,
												const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice)
:ARM_CRACalculator(mktDataManager),
itsLocalResetFreq(0),
itsLocalModelType(LocalDownUp), //DownUp for spread case !
itsLocalPortfolio(ARM_StdPortfolioPtr(NULL)),
itsIsExerciseProbas(false),
itsExerciseProbaOffset(0),
itsExerciseProbas(std::vector<double>(1,0.0))
{
	SetDegeneratedCalculator(CRAFromCCSO);

	Set(optPortfolio, 
		0.0, 0.0, 0.0, 0.0, //MRS and Beta = 0
		NULL,				//calibSecPFParams
		1,					//nbSteps
		1,					//flagToGenerateOSWATM
		mdmKeys, 
		productsToPrice, 
		-1,					//recalib Mrs
		-1,					//recalib Beta
		K_DIAG, 
		true);

	CheckCRAInputs();

	CheckData();

	GenerateProductDescription((ARM_OptionPortfolio *) optPortfolio);

	ARM_StringVector pricedColumns = PricedColumnNames();

	ARM_CstManagerPtr cstManagerPtr = CreateCstManager();

	Initialize();

	CreateAndSetDealDescriptionAndTimeIt("LOCALCAPFLOOR", pricedColumns, cstManagerPtr);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRALocalCalculator 
///	Routine: Summit Constructor from Option Portfolio 
///	Returns: void
///	Action : builds the object
/// WARNING: this constructor is called ONLY from CRA Spread
///          to degenerate the CCSO into CRA or Swaption
///          then no need to initialise the SFRM and local models
/////////////////////////////////////////////////////////////////
ARM_CRALocalCalculator::ARM_CRALocalCalculator( const ARM_Date& asOfDate,
												ARM_Security* security)
:ARM_CRACalculator(asOfDate, security),
itsLocalResetFreq(0),
itsLocalModelType(LocalDownUp), //DownUp for spread case !
itsLocalPortfolio(ARM_StdPortfolioPtr(NULL)),
itsIsExerciseProbas(false),
itsExerciseProbaOffset(0),
itsExerciseProbas(std::vector<double>(1,0.0))
{
	CheckCRAInputs();

	if (security->GetName() == ARM_OPTIONPORTFOLIO)
	{
		SetDegeneratedCalculator(CRAFromCCSO);

		ARM_OptionPortfolio* optionPf = (ARM_OptionPortfolio*) security;
		
		// Keep UP Barrier down
		// Force Down Barrier to -100 if 0
		optionPf->GetCorridorLeg()->GetDownBarriers()->SetCalcMethod(K_STEPUP_RIGHT);
		*(optionPf->GetCorridorLeg()->GetDownBarriers()) *= 100.0;
		ARM_Vector* BarDownValues = optionPf->GetCorridorLeg()->GetDownBarriers()->GetDiscreteValues();
    
		int bDownSz = BarDownValues->GetSize();
		for (int i = 0; i < bDownSz; i++)
		{
			if (fabs(BarDownValues->Elt(i)-0.0) < 1e-8)
				BarDownValues->Elt(i) = -10000.0;
		}

		// Keep UP Barrier
		optionPf->GetCorridorLeg()->GetUpBarriers()->SetCalcMethod(K_STEPUP_RIGHT);
		*(optionPf->GetCorridorLeg()->GetUpBarriers()) *= 100.0;

		// Funding spread, interpolation on Start dates
		ARM_ReferenceValue* theSpreads = optionPf->GetLiborLeg()->GetSpreads();

		ARM_Vector* resetDates = optionPf->GetLiborLeg()->GetResetDates();
		ARM_Vector* stDates =  (ARM_Vector *) optionPf->GetLiborLeg()->GetFlowStartDates()->Clone();

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
			ARM_Vector* resSpread = new ARM_Vector(resetDates->GetSize(), optionPf->GetLiborLeg()->GetSpread());
			ARM_ReferenceValue genSpreads(stDates, resSpread);	
			genSpreads.SetCalcMethod(K_LINEAR);
			optionPf->GetLiborLeg()->SetSpreads(&genSpreads);
			theSpreads = optionPf->GetLiborLeg()->GetSpreads();
		}
		*theSpreads *= 100.0;

		// Boosted Fix Rate
		optionPf->GetCorridorLeg()->GetSpreads()->SetCalcMethod(K_STEPUP_RIGHT);
		*(optionPf->GetCorridorLeg()->GetSpreads()) *= 100.0;

		GenerateProductDescription((ARM_OptionPortfolio*) optionPf);
	}
	else //it is a SWAPTION:
	{
		SetDegeneratedCalculator(SwaptionFromCCSO);

		GenerateProductDescription((ARM_Swaption*) security);
	}

	ARM_StringVector pricedColumns = PricedColumnNames();

	ARM_CstManagerPtr cstManagerPtr = CreateCstManager();

	CreateAndSetDealDescriptionAndTimeIt("LOCALCAPFLOOR", pricedColumns, cstManagerPtr);
}

////////////////////////////////////////////////////////////////////
///	Class  : ARM_CRALocalCalculator 
///	Routine: Simple constructor to transfer data to CRA constructor
///	Returns: void
///	Action : builds the object 
////////////////////////////////////////////////////////////////////
ARM_CRALocalCalculator::ARM_CRALocalCalculator(
					  const ARM_Currency& ccy,
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
					  bool isPortfolioNoticeDays,
					  int localModelType,
					  const ARM_StringVector& mdmKeys,
					  const ARM_MarketData_ManagerRep& mktDataManager,
				      const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice,
					  int reCalibMrs,
					  int reCalibBeta,
					  bool isStdCalib,
					  bool isExerciseProbas,
					  size_t exerciseProbaOffset,
					  ARM_Currency* fundccy,
					  const ARM_ReferenceValue& fundNominal )   							  
:ARM_CRACalculator(
ccy,
startDate,
endDate,
payRec,
notional,
callFreq,
callNotice,
callCal,
callFees,
fundFreq,
fundSpread,
fundDayCount,
cpnSpread,
cpnPayFreq,
cpnResetCal,
cpnPayCal,
boostedIndexType,
boostedFixRate,
boostedVarTerm,
boostedResetGap,
boostedResetTiming,
boostedDayCount,
boostedAdjRule,
boostedIntRule,
cpnBarrierDown,
cpnBarrierUp,
cpnResetFreq,
cpnResetTiming,
cpnResetGap,
refIndexType1,
refTerm1,
refDayCount1,
refCoeff1,
isPortfolioNoticeDays,
mdmKeys,
mktDataManager,
productsToPrice,
reCalibMrs,
reCalibBeta,
isStdCalib,
fundccy,
fundNominal),

itsLocalModelType(localModelType),
itsLocalResetFreq(0),
itsIsExerciseProbas(isExerciseProbas),
itsExerciseProbaOffset(exerciseProbaOffset),
itsExerciseProbas(std::vector<double>(1,0.0))
{
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRALocalCalculator
///	Routine: Copy constructor
///	Returns: a built object
///	Action : builds the object
/////////////////////////////////////////////////////////////////
ARM_CRALocalCalculator::ARM_CRALocalCalculator(const ARM_CRALocalCalculator& rhs)
					:	ARM_CRACalculator(rhs)
{
	CopyNoCleanUp(rhs);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRALocalCalculator
///	Routine: Destructor
///	Returns: void
///	Action : destroys the object
/////////////////////////////////////////////////////////////////
ARM_CRALocalCalculator::~ARM_CRALocalCalculator()
{
	CleanUp();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRALocalCalculator
///	Routine: assignment operator
///	Returns: 
///	Action : 
/////////////////////////////////////////////////////////////////
ARM_CRALocalCalculator& ARM_CRALocalCalculator::operator=(const ARM_CRALocalCalculator& rhs)
{
	if( this != & rhs )
	{
		ARM_CRACalculator::operator=(rhs);
        CleanUp();
        CopyNoCleanUp(rhs);
	}
	return *this;
}

////////////////////////////////////////////////////
///	Class   : ARM_CRALocalCalculator
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
///           Call copy constructor
////////////////////////////////////////////////////
ARM_Object* ARM_CRALocalCalculator::Clone() const
{
	return(new ARM_CRALocalCalculator(*this));
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRALocalCalculator
///	Routine: UpdateCalibration
///	Returns: void
///	Action : update the calibration datas w.r.t. MarketDataManager
///          The context is an hedge ratios computation
/////////////////////////////////////////////////////////////////
void ARM_CRALocalCalculator::UpdateCalibration(bool isUpdateStrike)
{
	ARM_CRACalculator::UpdateCalibration();

	ARM_BSSmiledModel* capModel = dynamic_cast< ARM_BSSmiledModel* >( GetMktDataManager()->GetData(GetKeys()[CfModelKey]) );
	/// Recompute market prices of caplets and swaptions		
	int i = 0;

	size_t size = itsLocalPortfolio->size();
	for (i = 0; i < size; ++i)
	{
		itsLocalPortfolio->GetAsset(i)->SetModelVariable(NULL);
		itsLocalPortfolio->GetAsset(i)->SetModel(capModel);

		double mktPrice = itsLocalPortfolio->GetAsset(i)->ComputePrice();
		
		itsLocalPortfolio->SetWeight(1.0, i);
		itsLocalPortfolio->SetPrecision(0.000001, i);
		itsLocalPortfolio->SetPrice(mktPrice, i);
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRALocalCalculator
///	Routine: ComputeExerciseDatesAndSetFees
///	Returns: void
///	Action : use this function instead of the CRA's one if 
///          called from a CCSO
/////////////////////////////////////////////////////////////////
std::vector<double>& ARM_CRALocalCalculator::ComputeExerciseDatesAndSetFees() const
{
	try
	{
		if ( GetDegeneratedCalculator() == No )
		   return ARM_CRACalculator::ComputeExerciseDatesAndSetFees();
		else
		{
			std::vector<double>& NewResetDates = NULL;

			if (itsOptionPortfolio)
			{
				ARM_CorridorLeg* corridor = (ARM_CorridorLeg*)itsOptionPortfolio->GetCorridorLeg();
				ARM_Vector* StartDates = corridor->GetFlowStartDates();
				ARM_Vector* ResetDates = corridor->GetResetDates();
				int legSize = ResetDates->GetSize();

				NewResetDates = To_pARM_GP_Vector(ResetDates);

				ARM_Vector* NoticeDates = itsOptionPortfolio->GetExerciseStyle()->GetExerciseStartDates();
				ARM_Vector* ExpiryDates = itsOptionPortfolio->GetExerciseStyle()->GetExerciseEndDates();
				int noticeSize = NoticeDates->GetSize();

				// be careful when notice frequency < payment frequency
				int ratio = 1;
				if ( (NoticeDates->size() > 1) && (StartDates->size() > 1) )
					ratio = floor(CountYears(corridor->GetDayCount(), NoticeDates->Elt(0), NoticeDates->Elt(1)) / CountYears(corridor->GetDayCount(), StartDates->Elt(0), StartDates->Elt(1)) + 0.5);

				int lastIdx = 0;
				for (int i = 0; i < legSize; i++)
				{
					for (int j = 0; j < noticeSize; j++)
					{
						if ( ARM_Date(ExpiryDates->Elt(j)) == ARM_Date(StartDates->Elt(i)) )
						{
							NewResetDates->Elt(i) = NoticeDates->Elt(j);
							break;
						}
					}
					// if not found, then create an intermediate date !
					if ( j == noticeSize )
					{
						if ( ratio < 1 )
						{
							char msg[200];
							sprintf(msg, "CraLocalCalculator::ComputeExerciseDatesAndSetFees : cannot generate exercise date for corridor's StartDate[%d]", i);
							throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg);
						}

						if ( ratio >= 1 )
						{
							for (int k = 0; (k < ratio-1) && (i+k < legSize); k++)
							{
								NewResetDates->Elt(i+k) = -1.0;
							}
						}
					}
				}

				itsCallFees = (*(itsOptionPortfolio->GetStrike()));
			}
			else
			{
			   throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
						"ARM_CRALocalCalculator::ComputeExerciseDatesAndSetFees : The option portfolio should not be null.");
			}

			return NewResetDates;
		}
	}

	catch(Exception& x)
	{
	    x.DebugPrint();
		throw x;
	}
	
	catch (...)
	{
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in CraLocalCalculator::ComputeExerciseDatesAndSetFees" );
	}
}


/////////////////////////////////////////////////////////////////////////////
///	Class  : ARM_CRALocalCalculator
///	Routine: DatesStructure
///	Returns: ARM_DateStripVector
///	Action : create the list of all event dates of the CRA. 
///			 The DateStripCombiner merges event dates of each
///          legs
/////////////////////////////////////////////////////////////////
ARM_DateStripCombiner ARM_CRALocalCalculator::DatesStructure() const
{
	try
	{
        if (itsExerDateStrip.IsNull())
        {
           return(ARM_CRACalculator::DatesStructure());
        }
        else
        {
           ARM_DateStripVector SchedVect(1, NULL);

           SchedVect[CALL_CRA_SCHED] = &(*itsExerDateStrip);

           return(SchedVect);
        }
	}

	catch(Exception& x)
	{
	    x.DebugPrint();

		throw x;
	}

	catch (...)
	{
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in CraLocalCalculator::DatesStructure" );
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRALocalCalculator
///	Routine: GetDescSize
///	Returns: size_t
///	Action : compute number of columns of deal description
/////////////////////////////////////////////////////////////////
size_t ARM_CRALocalCalculator::GetDescSize() const
{
	return itsLocalModelType == LocalDownUp
		? sizeof(LocalCapFloorCRAColNamesTable)/sizeof(LocalCapFloorCRAColNamesTable[0])
		: sizeof(LocalFloorCRAColNamesTable)/sizeof(LocalFloorCRAColNamesTable[0]);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRALocalCalculator
///	Routine: MiddleRows
///	Returns: ARM_RowInfo
///	Action : create a row of a deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_CRALocalCalculator::MiddleRows(size_t eventIdx, 
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
		size_t descSize = GetDescSize(); 
	
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

		//Weights for exotic swap
		double wCorridor = 1.0;
		double wFunding  = 1.0;

		//NOTICE DATE
		double noticeDate = (*(datesStructure.GetDateStrip(CALL_CRA_SCHED)->GetResetDates()))[eventIdx];

		CC_Ostringstream noticeDateDesc;
		noticeDateDesc << CC_NS(std,fixed) << noticeDate;
		rowDescVec[localCfResetDate] = noticeDateDesc.str();
		rowTypeVec[localCfResetDate] = ARM_DATE_TYPE;

		//START DATE / CALL DATE
		double startDate = (*(datesStructure.GetDateStrip(CALL_CRA_SCHED)->GetFlowStartDates()))[eventIdx] ;
		CC_Ostringstream startDateDesc;
		startDateDesc << CC_NS(std,fixed) << startDate;
		rowDescVec[localCfStartDate] = startDateDesc.str();
		rowTypeVec[localCfStartDate] = ARM_DATE_TYPE;
		
		//MATURITY DATE
		CC_Ostringstream maturityDateDesc;
		maturityDateDesc << CC_NS(std,fixed) << maturityDate;
		rowDescVec[localCfMaturityDate] = maturityDateDesc.str();
		rowTypeVec[localCfMaturityDate] = ARM_DATE_TYPE;

		//FEES
		double fees = const_cast< ARM_ReferenceValue& >(itsCallFees).Interpolate(startDate);
		CC_Ostringstream feesDesc;
		feesDesc << CC_NS(std,fixed) << fees;
		rowDescVec[localCfFees] = feesDesc.str();
		rowTypeVec[localCfFees] = ARM_DOUBLE;

		//FUNDING LEG
		CC_Ostringstream fundingLegDesc;
		fundingLegDesc << "SWAP(" << ccy << ", StartDate[i], MaturityDate[i], 0," << payRec << ", ";
		fundingLegDesc << fundFreq << ", " << fundDayCount << ", " << fundFreq << ", " << fundDayCount << ", ";
		fundingLegDesc << "FundSpreadCst, NotionalCst" << ")";
		rowDescVec[localCfFundingLeg] = fundingLegDesc.str();
		rowTypeVec[localCfFundingLeg] = ARM_STRING;


		//Local Cap Floor Deal Description
		if (itsLocalModelType == LocalDownUp)
		{
			//CORRIDOR LEG
			CC_Ostringstream corridorLegDesc;
			corridorLegDesc << "CORRIDOR(" << "LOCALCAPFLOOR" << ",";			// Model
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
			corridorLegDesc << cpnResetGap << ",";					// Pay Gap
			corridorLegDesc << cpnIntRule << ",";					// Int Rule
			corridorLegDesc << "CpnSpreadCst" << ",";				// Spread
			corridorLegDesc <<  cpnResetFreq << ",";				// Fixing Freq
			corridorLegDesc << cpnResetTiming << ",";				// Fixing Timing
			corridorLegDesc << cpnIndexType1 << ",";				// Fixing Index Type 1
			corridorLegDesc << cpnTerm1 << ",";						// Fixing Index Term 1
			corridorLegDesc << refCoeff << ",";				// Coeff 1
			corridorLegDesc << "BarrierDownCst" << ",";				// Barrier Down
			corridorLegDesc << "BarrierUpCst" << ",";				// Barrier Up
			corridorLegDesc << "NotionalCst" << ")";				// Notional
		
			rowDescVec[localCfCorridorLeg] = corridorLegDesc.str();
			rowTypeVec[localCfCorridorLeg] = ARM_STRING;

			//EXOTIC SWAP
			CC_Ostringstream ExoticSwapDesc;
			if (itsOptionPortfolio)
			{
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
			}

			ExoticSwapDesc << wCorridor << "*CorridorLeg[i]+" << wFunding << "*FundingLeg[i]";
			rowDescVec[localCfExoticSwap] = ExoticSwapDesc.str();
			rowTypeVec[localCfExoticSwap] = ARM_STRING;   

			//OPTION  
			CC_Ostringstream optionDesc;
			if(eventIdx == callEventSize-1)
			{
				optionDesc << "MAX(ExoticSwap[i]-Fees[i],0)";
			}
			else
			{
				optionDesc << "Exercise(0,ExoticSwap[i]-Fees[i],Option[i+1])";
			}
			rowDescVec[localCfOption] = optionDesc.str();
			rowTypeVec[localCfOption] = ARM_STRING;

			//BERMUDA
			CC_Ostringstream bermudaDesc;
			if(eventIdx == 0)
			{
				int payRec = itsOptionPortfolio ? itsOptionPortfolio->GetPorS() : itsSwaption->GetOptionType();
				if ( payRec == K_PAY)
				{
					bermudaDesc << "-Option[i]";
				}
				else
				{
					bermudaDesc << "Option[i]";
				}
				rowDescVec[localCfBermuda] = bermudaDesc.str();
				rowTypeVec[localCfBermuda] = ARM_STRING;
			}
			else
			{
				bermudaDesc << "0";
				rowDescVec[localCfBermuda] = bermudaDesc.str();
				rowTypeVec[localCfBermuda] = ARM_DOUBLE;
			}

			CC_Ostringstream corridorDesc;
			CC_Ostringstream fundingDesc;
			if(eventIdx == 0)
			{
				//CORRIDOR 
				corridorDesc << "CorridorLeg[i]";
				rowDescVec[localCfCorridor] = corridorDesc.str();
				rowTypeVec[localCfCorridor] = ARM_STRING;

				//FUNDING 
				fundingDesc << "FundingLeg[i]";
				rowDescVec[localCfFunding] = fundingDesc.str();
				rowTypeVec[localCfFunding] = ARM_STRING;
			}
			else
			{
				corridorDesc << "0";
				rowDescVec[localCfCorridor] = corridorDesc.str();
				rowTypeVec[localCfCorridor] = ARM_DOUBLE;
		
				fundingDesc << "0";
				rowDescVec[localCfFunding] = fundingDesc.str();
				rowTypeVec[localCfFunding] = ARM_DOUBLE;
			}

			//FWD CORRIDOR / FUNDING
			CC_Ostringstream fwdCorridorDesc;
			CC_Ostringstream fwdFundingDesc;
			if (eventIdx == fwdStart)
			{
				fwdCorridorDesc << "CorridorLeg[i]";
				rowDescVec[localCfFwdCorridor] = fwdCorridorDesc.str();
				rowTypeVec[localCfFwdCorridor] = ARM_STRING;
			
				fwdFundingDesc << "FundingLeg[i]";
				rowDescVec[localCfFwdFunding] = fwdFundingDesc.str();
				rowTypeVec[localCfFwdFunding] = ARM_STRING;
			}
			else
			{
				fwdCorridorDesc << "0";
				rowDescVec[localCfFwdCorridor] = fwdCorridorDesc.str();
				rowTypeVec[localCfFwdCorridor] = ARM_DOUBLE;

				CC_Ostringstream fwdFundingDesc;
				fwdFundingDesc << "0";
				rowDescVec[localCfFwdFunding] = fwdFundingDesc.str();
				rowTypeVec[localCfFwdFunding] = ARM_DOUBLE;
			}
		}
		//Local Floor / Local Cap Deal Description
		else
		{
			//CORRIDOR DOWN
			char* localDownModel; 
			if (itsLocalModelType == LocalDown)
				localDownModel = "LOCALFLOOR"; 
			else 
				localDownModel = ccy;
			CC_Ostringstream corridorDownDesc;
			corridorDownDesc << "CORRIDOR(" << localDownModel << ",";	// Model
			corridorDownDesc << "StartDate[i]" << ",";					// StartDate
			corridorDownDesc << "MaturityDate[i]" << ",";				// EndDate
			corridorDownDesc << payRec << ",";							// PayRec
			corridorDownDesc << boostedIndexType << ",";				// PayIndexType
			corridorDownDesc << "BoostedFixRateCst" << ",";				// Fix Value
			corridorDownDesc << 1 << ",";								// PayIndexMult Value
			corridorDownDesc << cpnPayFreq << ",";						// Payment Freq
			corridorDownDesc << boostedVarTerm << ",";					// Pay Index Term
			corridorDownDesc << boostedResetTiming << ",";				// Pay Index Timing
			corridorDownDesc << boostedPayDayCount << ",";				// Pay Day Count
			corridorDownDesc << cpnResetGap << ",";						// Pay Gap
			corridorDownDesc << cpnIntRule << ",";						// Int Rule
			corridorDownDesc << "CpnSpreadCst" << ",";					// Spread
			corridorDownDesc <<  cpnResetFreq << ",";					// Fixing Freq
			corridorDownDesc << cpnResetTiming << ",";					// Fixing Timing
			corridorDownDesc << cpnIndexType1 << ",";					// Fixing Index Type 1
			corridorDownDesc << cpnTerm1 << ",";						// Fixing Index Term 1
			corridorDownDesc << refCoeff	 << ",";					// Coeff 1
			corridorDownDesc << 0 << ",";								// Barrier Down set to zero
			corridorDownDesc << "BarrierDownCst" << ",";				// Barrier Up set to barrier down
			corridorDownDesc << "NotionalCst" << ")";					// Notional
		
			rowDescVec[localFCorridorDown] = corridorDownDesc.str();
			rowTypeVec[localFCorridorDown] = ARM_STRING;

			//CORRIDOR UP
			char* localUpModel; 
			if (itsLocalModelType == LocalDown)
				localUpModel = ccy; 
			else 
				localUpModel = "LOCALCAP";
			CC_Ostringstream corridorUpDesc;
			corridorUpDesc << "CORRIDOR(" << localUpModel << ",";	// Model
			corridorUpDesc << "StartDate[i]" << ",";				// StartDate
			corridorUpDesc << "MaturityDate[i]" << ",";				// EndDate
			corridorUpDesc << payRec << ",";						// PayRec
			corridorUpDesc << boostedIndexType << ",";				// PayIndexType
			corridorUpDesc << "BoostedFixRateCst" << ",";			// Fix Value
			corridorUpDesc << 1				 << ",";				// PayIndexMult Value
			corridorUpDesc << cpnPayFreq << ",";					// Payment Freq
			corridorUpDesc << boostedVarTerm << ",";				// Pay Index Term
			corridorUpDesc << boostedResetTiming << ",";			// Pay Index Timing
			corridorUpDesc << boostedPayDayCount << ",";			// Pay Day Count
			corridorUpDesc << cpnResetGap << ",";						// Pay Gap
			corridorUpDesc << cpnIntRule << ",";					// Int Rule
			corridorUpDesc << "CpnSpreadCst" << ",";				// Spread
			corridorUpDesc <<  cpnResetFreq << ",";					// Fixing Freq
			corridorUpDesc << cpnResetTiming << ",";				// Fixing Timing
			corridorUpDesc << cpnIndexType1 << ",";					// Fixing Index Type 1
			corridorUpDesc << cpnTerm1 << ",";						// Fixing Index Term 1
			corridorUpDesc << refCoeff			<< ",";				// Coeff 1
			corridorUpDesc << 0 << ",";								// Barrier Down set to Zero
			corridorUpDesc << "BarrierUpCst" << ",";				// Barrier Up
			corridorUpDesc << "NotionalCst" << ")";					// Notional
		
			rowDescVec[localFCorridorUp] = corridorUpDesc.str();
			rowTypeVec[localFCorridorUp] = ARM_STRING;
			
			//CORRIDOR
			CC_Ostringstream corridorLegDesc;
			corridorLegDesc << "CorridorUp[i]-CorridorDown[i]";
			rowDescVec[localFCorridorLeg] = corridorLegDesc.str();
			rowTypeVec[localFCorridorLeg] = ARM_STRING;

			//EXOTIC SWAP
			CC_Ostringstream ExoticSwapDesc;
			ExoticSwapDesc << wCorridor << "*CorridorLeg[i]+" << wFunding << "*FundingLeg[i]";
			rowDescVec[localFExoticSwap] = ExoticSwapDesc.str();
			rowTypeVec[localFExoticSwap] = ARM_STRING;   
		
			//OPTION  
			CC_Ostringstream optionDesc;
			if(eventIdx == callEventSize-1)
			{
				optionDesc << "MAX(ExoticSwap[i]-Fees[i],0)";
			}
			else
			{
				optionDesc << "Exercise(0,ExoticSwap[i]-Fees[i],Option[i+1])";
			}
			rowDescVec[localFOption] = optionDesc.str();
			rowTypeVec[localFOption] = ARM_STRING;

			//BERMUDA
			CC_Ostringstream bermudaDesc;
			if(eventIdx == 0)
			{
				if (itsOptionPortfolio->GetPorS() == K_PAY)
				{
					bermudaDesc << "-Option[i]";
				}
				else
				{
					bermudaDesc << "Option[i]";
				}
				rowDescVec[localFBermuda] = bermudaDesc.str();
				rowTypeVec[localFBermuda] = ARM_STRING;
			}
			else
			{
				bermudaDesc << "0";
				rowDescVec[localFBermuda] = bermudaDesc.str();
				rowTypeVec[localFBermuda] = ARM_DOUBLE;
			}
			
			CC_Ostringstream corridorDesc;	
			CC_Ostringstream fundingDesc;
			if(eventIdx == 0)
			{
				//CORRIDOR 
				corridorDesc << "CorridorLeg[i]";
				rowDescVec[localFCorridor] = corridorDesc.str();
				rowTypeVec[localFCorridor] = ARM_STRING;

				//FUNDING 
				fundingDesc << "SWAP(" << ccy << ", StartDate[i], MaturityDate[i], 0," << payRec << ", ";
				fundingDesc << fundFreq << ", " << fundDayCount << ", " << fundFreq << ", " << fundDayCount << ", ";
				fundingDesc << "FundSpreadCst, NotionalCst" << ")";
				rowDescVec[localFFunding] = fundingDesc.str();
				rowTypeVec[localFFunding] = ARM_STRING;
			}
			else
			{
				corridorDesc << "0";
				rowDescVec[localFCorridor] = corridorDesc.str();
				rowTypeVec[localFCorridor] = ARM_DOUBLE;
		
				fundingDesc << "0";
				rowDescVec[localFFunding] = fundingDesc.str();
				rowTypeVec[localFFunding] = ARM_DOUBLE;
			}
			
			//FWD CORRIDOR / FUNDING
			CC_Ostringstream fwdCorridorDownDesc;
			CC_Ostringstream fwdCorridorUpDesc;
			CC_Ostringstream fwdCorridorDesc;
			CC_Ostringstream fwdFundingDesc;
			if (eventIdx == fwdStart)
			{
				fwdCorridorDownDesc << "CORRIDOR(" << localDownModel << ",";	// Model
				fwdCorridorDownDesc << "StartDate[i]" << ",";					// StartDate
				fwdCorridorDownDesc << "MaturityDate[i]" << ",";				// EndDate
				fwdCorridorDownDesc << payRec << ",";							// PayRec
				fwdCorridorDownDesc << boostedIndexType << ",";					// PayIndexType
				fwdCorridorDownDesc << "BoostedFixRateCst" << ",";				// Fix Value
				fwdCorridorDownDesc << 1				 << ",";				// PayIndexMult Value
				fwdCorridorDownDesc << cpnPayFreq << ",";						// Payment Freq
				fwdCorridorDownDesc << boostedVarTerm << ",";					// Pay Index Term
				fwdCorridorDownDesc << boostedResetTiming << ",";				// Pay Index Timing
				fwdCorridorDownDesc << boostedPayDayCount << ",";				// Pay Day Count
				fwdCorridorDownDesc << cpnResetGap << ",";						// Pay Gap
				fwdCorridorDownDesc << cpnIntRule << ",";						// Int Rule
				fwdCorridorDownDesc << "CpnSpreadCst" << ",";					// Spread
				fwdCorridorDownDesc <<  cpnResetFreq << ",";					// Fixing Freq
				fwdCorridorDownDesc << cpnResetTiming << ",";					// Fixing Timing
				fwdCorridorDownDesc << cpnIndexType1 << ",";					// Fixing Index Type 1
				fwdCorridorDownDesc << cpnTerm1 << ",";							// Fixing Index Term 1
				fwdCorridorDownDesc << refCoeff		 << ",";					// Coeff 1
				fwdCorridorDownDesc << 0 << ",";								// Barrier Down set to zero
				fwdCorridorDownDesc << "BarrierDownCst" << ",";					// Barrier Up set to barrier down
				fwdCorridorDownDesc << "NotionalCst" << ")";					// Notional
			
				rowDescVec[localFFwdCorridorDown] = fwdCorridorDownDesc.str();
				rowTypeVec[localFFwdCorridorDown] = ARM_STRING;

				//FWD CORRIDOR UP
				fwdCorridorUpDesc << "CORRIDOR(" << localUpModel << ",";	// Model
				fwdCorridorUpDesc << "StartDate[i]" << ",";					// StartDate
				fwdCorridorUpDesc << "MaturityDate[i]" << ",";				// EndDate
				fwdCorridorUpDesc << payRec << ",";							// PayRec
				fwdCorridorUpDesc << boostedIndexType << ",";				// PayIndexType
				fwdCorridorUpDesc << "BoostedFixRateCst" << ",";			// Fix Value
				fwdCorridorUpDesc << 1				 << ",";				// PayIndexMult Value
				fwdCorridorUpDesc << cpnPayFreq << ",";						// Payment Freq
				fwdCorridorUpDesc << boostedVarTerm << ",";					// Pay Index Term
				fwdCorridorUpDesc << boostedResetTiming << ",";				// Pay Index Timing
				fwdCorridorUpDesc << boostedPayDayCount << ",";				// Pay Day Count
				fwdCorridorUpDesc << cpnResetGap << ",";					// Pay Gap
				fwdCorridorUpDesc << cpnIntRule << ",";						// Int Rule
				fwdCorridorUpDesc << "CpnSpreadCst" << ",";					// Spread
				fwdCorridorUpDesc <<  cpnResetFreq << ",";					// Fixing Freq
				fwdCorridorUpDesc << cpnResetTiming << ",";					// Fixing Timing
				fwdCorridorUpDesc << cpnIndexType1 << ",";					// Fixing Index Type 1
				fwdCorridorUpDesc << cpnTerm1 << ",";						// Fixing Index Term 1
				fwdCorridorUpDesc << refCoeff	 << ",";					// Coeff 1
				fwdCorridorUpDesc << 0 << ",";								// Barrier Down set to Zero
				fwdCorridorUpDesc << "BarrierUpCst" << ",";					// Barrier Up
				fwdCorridorUpDesc << "NotionalCst" << ")";					// Notional
			
				rowDescVec[localFFwdCorridorUp] = fwdCorridorUpDesc.str();
				rowTypeVec[localFFwdCorridorUp] = ARM_STRING;	
				
				//FWD CORRIDOR
				fwdCorridorDesc << "FwdCorridorUp[i]-FwdCorridorDown[i]";
				rowDescVec[localFFwdCorridor] = fwdCorridorDesc.str();
				rowTypeVec[localFFwdCorridor] = ARM_STRING;

				//FWD FUNDING
				fwdFundingDesc << "SWAP(" << ccy << ", StartDate[i], MaturityDate[i], 0," << payRec << ", ";
				fwdFundingDesc << fundFreq << ", " << fundDayCount << ", " << fundFreq << ", " << fundDayCount << ", ";
				fwdFundingDesc << "FundSpreadCst, NotionalCst" << ")";
				rowDescVec[localFFwdFunding] = fwdFundingDesc.str();
				rowTypeVec[localFFwdFunding] = ARM_STRING;
			}
			else
			{
				fwdCorridorUpDesc << "0";
				rowDescVec[localFFwdCorridorUp] = fwdCorridorUpDesc.str();
				rowTypeVec[localFFwdCorridorUp] = ARM_DOUBLE;

				fwdCorridorDownDesc << "0";
				rowDescVec[localFFwdCorridorDown] = fwdCorridorDownDesc.str();
				rowTypeVec[localFFwdCorridorDown] = ARM_DOUBLE;

				fwdCorridorDesc << "0";
				rowDescVec[localFFwdCorridor] = fwdCorridorDesc.str();
				rowTypeVec[localFFwdCorridor] = ARM_DOUBLE;

				CC_Ostringstream fwdFundingDesc;
				fwdFundingDesc << "0";
				rowDescVec[localFFwdFunding] = fwdFundingDesc.str();
				rowTypeVec[localFFwdFunding] = ARM_DOUBLE;
			}
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
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in CraLocalCalculator::MiddleRows" );
	}
}

////////////////////////////////////////////////////////////////////////////
///	Class  : ARM_CRALocalCalculator
///	Routine: PricedColumnNames
///	Returns: ARM_StringVector
///	Action : create the priced column names of the deal description
//////////////////////////////////////////////////////////////////////////////
ARM_StringVector ARM_CRALocalCalculator::PricedColumnNames() const
{
	ARM_StringVector pricedColumns;

	int colSize = 0;
	int i = 0;

	//LocalCapFloor
	if (itsLocalModelType == LocalDownUp)
	{
		// effective nb of columns to price :
		colSize = MIN(localCfNbProductsToPrice, itsProductsToPrice.size());
	}
	//LocalFloor
	else
	{	
		colSize = MIN(localFNbProductsToPrice, itsProductsToPrice.size());
	}

	//LocalCapFloor
	if (itsLocalModelType == LocalDownUp)
	{
		for (i=0; i<colSize; i++)
		{
			if (itsProductsToPrice[i])
				pricedColumns.push_back(LocalCapFloorCRAColNamesTable[LocalCapFloorCRAProductToPriceColumns[i]]);
		}

		if(itsIsExerciseProbas)
		{
			for(size_t i=0;i<LocalCapFloorCRAExerciseProbaNb;++i)
				pricedColumns.push_back(LocalCapFloorCRAExerciseProbaNamesTable[i]);
		}
	}
	//LocalFloor
	else
	{	
		for (i=0; i<colSize; i++)
		{
			if (itsProductsToPrice[i])
				pricedColumns.push_back(LocalFloorCRAColNamesTable[LocalFloorCRAProductToPriceColumns[i]]);
		}
	}

	return pricedColumns;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRALocalCalculator
///	Routine: ColumnNames
///	Returns: ARM_RowInfo
///	Action : create the column names of the deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_CRALocalCalculator::ColumnNames() const
{
	// Number of Columns
	size_t colNamesSize = GetDescSize();

    vector< string > colNamesVec(colNamesSize);
    vector< ARM_GP_VALUE_TYPE > colTypeVec(colNamesSize, ARM_STRING); 

	//LocalCapFloor
	if (itsLocalModelType == LocalDownUp)
	{
		for (size_t i = 0; i < colNamesSize; ++i)
		{
			colNamesVec[i] = LocalCapFloorCRAColNamesTable[i];
		}
	}
	//LocalFloor
	else
	{
		for (size_t i = 0; i < colNamesSize; ++i)
		{
			colNamesVec[i] = LocalFloorCRAColNamesTable[i];
		}
	}

    ARM_RowInfo rowInfo(colNamesVec,colTypeVec);

    return rowInfo;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRALocalCalculator
///	Routine: InitPriceableColumns
///	Returns: nothing
///	Action : Feed the MDM with input market objects
/////////////////////////////////////////////////////////////////
void ARM_CRALocalCalculator::InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const
{
    string zeroValue("0");

	if (itsLocalModelType == LocalDownUp)
	{
		rowDescVec[localCfBermuda] = zeroValue;
		rowTypeVec[localCfBermuda] = ARM_DOUBLE;

		rowDescVec[localCfFundingLeg] = zeroValue;
		rowTypeVec[localCfFundingLeg] = ARM_DOUBLE;

		rowDescVec[localCfCorridorLeg] = zeroValue;
		rowTypeVec[localCfCorridorLeg] = ARM_DOUBLE;

		rowDescVec[localCfExoticSwap] = zeroValue;
		rowTypeVec[localCfExoticSwap] = ARM_DOUBLE;

		rowDescVec[localCfExerSwapRate] = zeroValue;
		rowTypeVec[localCfExerSwapRate] = ARM_DOUBLE;
		
		rowDescVec[localCfFrontier] = zeroValue;
		rowTypeVec[localCfFrontier] = ARM_DOUBLE;

		for(size_t i=0;i<LocalCapFloorCRAExerciseProbaNb;++i)
		{
			rowDescVec[localCfProba1+2*i] = zeroValue;
			rowTypeVec[localCfProba1+2*i] = ARM_DOUBLE;
		}
	}
	else
	{
		rowDescVec[localFBermuda] = zeroValue;
		rowTypeVec[localFBermuda] = ARM_DOUBLE;

		rowDescVec[localFFundingLeg] = zeroValue;
		rowTypeVec[localFFundingLeg] = ARM_DOUBLE;

		rowDescVec[localFCorridorLeg] = zeroValue;
		rowTypeVec[localFCorridorLeg] = ARM_DOUBLE;

		rowDescVec[localFExoticSwap] = zeroValue;
		rowTypeVec[localFExoticSwap] = ARM_DOUBLE;
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRALocalCalculator
///	Routine: CreateAndSetModel
///	Returns: void
///	Action : create the encapsulated model
/////////////////////////////////////////////////////////////////
void ARM_CRALocalCalculator::CreateAndSetModel()
{
	try
	{
		ARM_CRACalculator::CreateAndSetModel();
		
		CreateAndSetLocalModel();
	}

	catch (Exception& x)
	{
	    x.DebugPrint();

		throw x;
	}

	catch (...)
	{
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in CRALocalCalculator::CreateAndSetModel" );
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRALocalCalculator
///	Routine: CreateAndSetLocalModel
///	Returns: void
///	Action : create the encapsulated local model
/////////////////////////////////////////////////////////////////
void ARM_CRALocalCalculator::CreateAndSetLocalModel()
{
	try
	{
		string localModelName;
		if (itsLocalModelType == LocalDown)
			localModelName = "LOCALFLOOR";
		else if (itsLocalModelType == LocalUp)
			localModelName = "LOCALCAP";
		else 
			localModelName = "LOCALCAPFLOOR";

		string ccy = string(itsCcy.GetCcyName());
		ARM_InterpolType type = ARM_InterpolationType::linear_column_extrapoleCst;

		//FORWARD ADJUSTMENT
		vector<double> fwdAdjCall = vector<double>(1);
		fwdAdjCall[0] = 0.0;
		vector<double> fwdAdjFix = vector<double>(1);
		fwdAdjFix[0] = 0.0;
		vector<double> fwdAdjValue = vector<double>(1);
		fwdAdjValue[0] = 1.0;
		std::vector<double> vecFwdAdjCall(fwdAdjCall);
		std::vector<double> vecFwdAdjFix(fwdAdjFix);
		ARM_GP_T_Matrix<double> matFwdAdjValue(1, 1, fwdAdjValue);
		ARM_SurfaceWithInterpol* fwdAdjSurface = new ARM_SurfaceWithInterpol( vecFwdAdjCall, vecFwdAdjFix, matFwdAdjValue, type);
		vector<double> fwdAdjIndex = vector<double>(1);
		fwdAdjIndex[0] = 0.0;
		ARM_SurfacePtrVector fwdAdjSurfaceList;  
		ARM_SurfacePtr fwdAdjSurfacePtr =  ARM_SurfacePtr( static_cast<ARM_Surface*>(fwdAdjSurface) );
		fwdAdjSurfaceList.push_back(fwdAdjSurfacePtr);
		ARM_ModelParam* fwdAdjModelParam = new ARM_SurfaceListModelParam(   ARM_ModelParamType::ForwardAdjustment,
																			std::vector<double>(fwdAdjIndex),
																			fwdAdjSurfaceList);															   
				
		//VOLATILITY
		int volSize = 4;
		vector<double> volCall = vector<double>(1);
		vector<double> volFix = vector<double>(1);
		vector<double> volValue = vector<double>(1);
		volCall[0] = 0.0;
		volFix[0] = 0.0;
		volValue[0] = 0.2;
		std::vector<double> vecVolCall(volCall);
		std::vector<double> vecVolFix(volFix);
		ARM_GP_T_Matrix<double> matVolValue(1, 1, volValue);
		ARM_SurfaceWithInterpol* volSurface = new ARM_SurfaceWithInterpol( vecVolCall, vecVolFix, matVolValue, type);
		vector<double> volIndex = vector<double>(4);
		volIndex[0] = 0.0;
		volIndex[1] = 1.0;
		volIndex[2] = 2.0;
		volIndex[3] = 3.0;
		ARM_SurfacePtrVector volSurfaceList;  
		
		for (int i = 0; i < volSize ; i++)
		{
			ARM_SurfacePtr volSurfacePtr =  ARM_SurfacePtr( static_cast<ARM_Surface*>(volSurface->Clone()) );
			volSurfaceList.push_back(volSurfacePtr);
		}
	
		ARM_ModelParam* volModelParam = new ARM_SurfaceListModelParam(  ARM_ModelParamType::Volatility,
																		std::vector<double>(volIndex),
																		volSurfaceList);
		if (volSurface)
		{
			delete volSurface;
			volSurface = NULL;
		}

		//SHIFT
		vector<double> shiftCall = vector<double>(1);
		shiftCall[0] = 0.0;
		vector<double> shiftFix = vector<double>(1);
		shiftFix[0] = 0.0;
		vector<double> shiftValue = vector<double>(1);
		shiftValue[0] = 0.0;
		std::vector<double> vecShiftCall(shiftCall);
		std::vector<double> vecShiftFix(shiftFix);
		ARM_GP_T_Matrix<double> matShiftValue(1, 1, shiftValue);
		ARM_SurfaceWithInterpol* shiftSurface = new ARM_SurfaceWithInterpol( vecShiftCall, vecShiftFix, matShiftValue, type);
		vector<double> shiftIndex = vector<double>(1);
		shiftIndex[0] = 0.0;
		ARM_SurfacePtrVector shiftSurfaceList;  
		ARM_SurfacePtr shiftSurfacePtr =  ARM_SurfacePtr( static_cast<ARM_Surface*>(shiftSurface) );
		shiftSurfaceList.push_back(shiftSurfacePtr);
		ARM_ModelParam* shiftModelParam = new ARM_SurfaceListModelParam(  ARM_ModelParamType::Shift,
																		  std::vector<double>(shiftIndex),
																		  shiftSurfaceList);
													
		//SLN MODEL PARAMS
		CC_STL_VECTOR( ARM_ModelParam* ) modelParams (3);
		modelParams[0] = fwdAdjModelParam;
		modelParams[1] = volModelParam;
		modelParams[2] = shiftModelParam;
		ARM_Local_SLN_ModelParams* localModelParams = new ARM_Local_SLN_ModelParams(modelParams);
		delete fwdAdjModelParam;
		delete volModelParam;
		delete shiftModelParam;

		//LOCAL SLN MODEL
		ARM_AutoCleaner<ARM_Local_SLN_ModelParams> Hold(localModelParams);
			
		//MODEL NAME MAP
		CC_STL_VECTOR( ARM_PricingModelPtr ) models(2);
		ARM_PricingModelPtr sfrmModel = GetPricingModel();
		
		models[0] = sfrmModel; 	
		models[1] = ARM_PricingModelPtr( new ARM_Local_SLN_Model( *localModelParams ) );
			
		vector<string> names = vector<string>(2);
		names[0] = ccy;
		names[1] = localModelName;

		ARM_StringVectorVector otherModelNames(2);
		otherModelNames[1] = ARM_StringVector (1,ccy); 
			
		ARM_ModelNameMap* modelNameMap = new ARM_ModelNameMap( names, models, otherModelNames );
				
		//CORREL MATRIX
		ARM_GP_Matrix* correlMatrix = new ARM_GP_Matrix(2, 2);
		(*correlMatrix)(0,0) = 1.0;
		(*correlMatrix)(1,0) = 0.0;
		(*correlMatrix)(0,1) = 0.0;
		(*correlMatrix)(1,1) = 1.0;

		//MULTI ASSET MODEL
		ARM_PricingModelPtr refModel (new ARM_MultiAssetsModel ( modelNameMap, correlMatrix ) );
		delete correlMatrix;

		//NUMERAIRE
		ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::TerminalZc ) );
		refModel->SetNumeraire(numeraire);
			
		//NUM METHOD
		int schedulerType = ARM_SchedulerBase::ConstantVariance;
		int samplerType = ARM_SamplerBase::NormalCentred;
		int truncatorType = ARM_TruncatorBase::StandardDeviation;
		int reconnectorType = ARM_ReconnectorBase::Mean;
		int smootherType = ARM_SmootherBase::DoNothing;
		std::vector<double> schedulerDatas(3);
		schedulerDatas[0] = GetNbSteps();
		schedulerDatas[1] = treeParams[1]; 

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

		refModel->SetNumMethod(ARM_NumMethodPtr(tree));
				
		SetPricingModel(refModel);
	}

	catch (Exception& x)
	{
	    x.DebugPrint();

		throw x;
	}

	catch (...)
	{
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in CRALocalCalculator::CreateAndSetLocalModel" );
	}
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRALocalCalculator
///	Routine: CalibrateLocalModel
///	Returns: void
///	Action : create the encapsulated local model
/////////////////////////////////////////////////////////////////
void ARM_CRALocalCalculator::CalibrateLocalModel()
{
	try
	{
		//LOCAL MODEL NAMES
		string localModelName;
		if (itsLocalModelType == LocalDown)
			localModelName = "LOCALFLOOR";
		else if (itsLocalModelType == LocalUp)
			localModelName = "LOCALCAP";
		else 
			localModelName = "LOCALCAPFLOOR";

		string ccy = string(itsCcy.GetCcyName());

		//GET MULTI ASSET MODEL
		ARM_PricingModelPtr refModel			= GetPricingModel();
		ARM_MultiAssetsModel*  multiAssetsModel = dynamic_cast<ARM_MultiAssetsModel*>(&*refModel);
	
		//LOCAL MODEL / SFRM MODEL
		ARM_Local_Model* calibratedLocalModel	= dynamic_cast<ARM_Local_Model*> (&*(*multiAssetsModel->GetModelMap())[localModelName]->Model());
		ARM_PricingModel* sfrmModel				= dynamic_cast<ARM_PricingModel*> (&*(*multiAssetsModel->GetModelMap())[ccy]->Model());

		//MULTI ASSET CALIBRATION = LOCAL MODEL CALIBRATION
		ARM_StdPortfolio* localPortfolio = &(*GetLocalPortfolio());
		int portfolioSize = localPortfolio->size();
		std::vector<double> evalTimes (portfolioSize);
		double currReset = 0.0;
		ARM_Date currDate;
		int noticeSize = GetExerDateStrip()->GetResetDates()->size();
		if (portfolioSize == noticeSize)
		{
			for (int i = 0; i < portfolioSize; i++) 
			{
				currReset = GetExerDateStrip()->GetResetDates()->Elt(i) ;
				currDate = ARM_Date(currReset);
				evalTimes[i] = sfrmModel->GetTimeFromDate(currDate);
			}
		}

		//CALIBRATION
		calibratedLocalModel->CalibrateLocalModel(*localPortfolio, evalTimes);	
	}

	catch (Exception& x)
	{
	    x.DebugPrint();

		throw x;
	}

	catch (...)
	{
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in CRALocalCalculator::CalibrateLocalModel" );
	}
}


////////////////////////////////////////////////////////////////////////////
///	Class  : ARM_CRALocalCalculator
///	Routine: CreateDownCorridorPortfolio
///	Returns: void
///	Action : create a corridor [0, bDown] portfolio for local calibration
////////////////////////////////////////////////////////////////////////////
void ARM_CRALocalCalculator::CreateCorridorLocalPortfolio()
{
	//Different barriers
	ARM_ReferenceValue barrierDown;
	ARM_ReferenceValue barrierUp;

	if (itsLocalModelType == LocalDown)
	{
		barrierDown = ARM_ReferenceValue(0.0);
		barrierUp = itsCpnBarrierDown;
	}
	else if (itsLocalModelType == LocalUp)
	{
		barrierDown = ARM_ReferenceValue(0.0);
		barrierUp = itsCpnBarrierUp;
	}
	else
	{
		barrierDown = itsCpnBarrierDown;
		barrierUp = itsCpnBarrierUp;
	}

	barrierDown *= 100.0;
	//barrierDown.SetCalcMethod(K_STEPUP_LEFT);
	barrierUp	*= 100.0;
	//barrierUp.SetCalcMethod(K_STEPUP_LEFT);

	//Deal description schedule.
	ARM_DateStripCombiner datesStructure	= DatesStructure();
	size_t portfolioSize = datesStructure.GetDateStrip(0)->GetResetDates()->size();
	double start; 
	ARM_Date startDate;
	ARM_Date nextStartDate;

	ARM_CorridorLeg* corridor = (ARM_CorridorLeg*)itsOptionPortfolio->GetCorridorLeg();
	//Some potential memory leaks.
	list < ARM_Security* > localCorridorList;
	ARM_IRIndex* paymentIndex			= NULL; 
	ARM_IRIndex* refIndex				= NULL;

	//Useful datas.
	int refResetFreq				= GetLocalResetFreq(); //GetCpnResetFreq();	
	int boostedResetTiming			= GetBoostedResetTiming();
	int cpnPayFreq					= GetCpnPayFreq();
	int cpnResetTiming				= GetCpnResetTiming();
	
	int addMonths; 
	if(cpnPayFreq == K_ANNUAL)
	{
		addMonths = 12;
	}
	if(cpnPayFreq == K_SEMIANNUAL)
	{
		addMonths = 6;
	}
	if(cpnPayFreq == K_QUARTERLY)
	{
		addMonths = 3;
	}
	if(cpnPayFreq == K_MONTHLY)
	{
		addMonths = 1;
	}

	ARM_ReferenceValue cpnSpread	= GetCpnSpread();
		
	if (itsBoostedIndexType == Fixed)
	{
		paymentIndex = new ARM_IRIndex(GetCcy().GetCcyName(), itsBoostedDayCount);
		cpnSpread = ARM_ReferenceValue(itsBoostedFixRate);
	}	
	else
	{
		paymentIndex= (ARM_IRIndex *)corridor->GetPaymentIndex()->Clone();
		
		paymentIndex->SetPayTiming(K_ARREARS);
		cpnSpread = ARM_ReferenceValue(itsBoostedFixRate);
	}
	cpnSpread *= 100.0;

	string liborTypeName(string(GetCcy().GetCcyName()) == string("EUR") ? "EURIBOR" : "LIBOR");
	liborTypeName += itsRefTerm1;
	ARM_INDEX_TYPE refIndexType = static_cast< ARM_INDEX_TYPE > (ARM_ArgConv_IndexType.GetNumber(liborTypeName));

	refIndex= (ARM_IRIndex *)corridor->GetRefIndex()->Clone(); 	

	if (refResetFreq != K_DEF_FREQ)
		refIndex->SetResetFrequency(refResetFreq);

	ARM_Currency ccy = GetCcy();

	ARM_ReferenceValue notional = GetNotional();
	//notional.SetCalcMethod(K_STEPUP_RIGHT);

	int payRec = GetPayRec();

	//Weights, precision
	double* weights		= new double[portfolioSize];
	double* precision	= new double[portfolioSize];
	double* mktPrices	= new double[portfolioSize];

	//Short term BS Smiled Model
	ARM_BSModel* capModel = dynamic_cast< ARM_BSModel* >(GetMktDataManager()->GetData(GetKeys()[CfModelKey]) );
		
	double mktPrice = 0.0;

	ARM_StdPortfolio* port = new ARM_StdPortfolio(portfolioSize, weights, mktPrices);

	int i = 0;

	int RefResetTiming = refIndex->GetResetTiming();
	ARM_Vector dateBarrier;
	if(RefResetTiming != K_ARREARS)
	{
		dateBarrier = *(corridor->GetPaymentStartDates());
	}
	else
	{
		dateBarrier = *(corridor->GetPaymentDates());
	}

	for (i = 0; i < portfolioSize; i++)
	{
		start = (*(datesStructure.GetDateStrip(0)->GetFlowStartDates()))[i] ;
		startDate = ARM_Date(start);
		nextStartDate = ARM_Date(start);
		nextStartDate.AddMonths(addMonths);
		double interpolDate = dateBarrier[i];
        double barrierDownValue = barrierDown.CptReferenceValue(interpolDate);
		ARM_ReferenceValue constbarrierDown(barrierDownValue);

		double barrierUpValue = barrierUp.CptReferenceValue(interpolDate);
		ARM_ReferenceValue constbarrierUp(barrierUpValue);

		ARM_CorridorLeg* corridorLet = new ARM_CorridorLeg(startDate, 
														   nextStartDate,
														   payRec, 
														   paymentIndex,
														   cpnPayFreq,
														   &cpnSpread,
														   refIndex, 
														   refResetFreq,
														   boostedResetTiming, 
														   cpnResetTiming,
														   K_SHORTSTART, //TMP
														   &constbarrierDown, 
														   K_STD,
														   &constbarrierUp, 
														   K_STD,
														   &ccy,
														   K_DEF_FREQ, 
														   K_LINEAR,
														   K_DIGITALE,
											 			   1, //TMP: forced CdecompPricingFlag,
														   (char*)(GetCpnResetCal()).c_str(),
														   (char*)(GetCpnPayCal()).c_str());
		

		corridorLet->SetModelVariable(NULL);
		corridorLet->SetModel(capModel);
		
		mktPrice = (corridorLet->ComputePrice());
			
		port->SetAsset(corridorLet, i);
		port->SetPrecision(0.000001,i);
		port->SetPrice(mktPrice, i);
		port->SetWeight(1.0, i);
	}

	//Free memory
	delete refIndex;
	refIndex = NULL;
	delete paymentIndex;
	paymentIndex = NULL;
	delete [] weights;
	weights = NULL;
	delete precision;
	precision = NULL;
	delete [] mktPrices;
	mktPrices = NULL;
	
	SetLocalPortfolio(ARM_StdPortfolioPtr(port)); 
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRALocalCalculator
///	Routine: UpdateModel
///	Returns: void
///	Action : update the model datas w.r.t. MarketDataManager
/////////////////////////////////////////////////////////////////
void ARM_CRALocalCalculator::UpdateModel()
{
	//ARM_CRACalculator::CreateAndSetModel();
	ARM_ZeroCurve* cpnCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));
	string ccy = string(itsCcy.GetCcyName());

	ARM_MultiAssetsModel*  multiAssetsModel = NULL;
	ARM_PricingModel* sfrmModel = NULL;

	multiAssetsModel = dynamic_cast<ARM_MultiAssetsModel*>(&*(GetPricingModel()));
	
	sfrmModel	 = dynamic_cast<ARM_PricingModel*> (&*(*multiAssetsModel->GetModelMap())[ccy]->Model());

	if( !sfrmModel)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : SFRM model is not a good type for updating");
	
	sfrmModel->SetZeroCurve( CreateClonedPtr( cpnCurve ));
		
	/// Update yield curves
	GetPricingModel()->SetZeroCurve( CreateClonedPtr( cpnCurve ) );
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRALocalCalculator
///	Routine: GetSFRMModel
///	Returns: ARM_PricingModel*
///	Action : get sfrm model from multi asset model
/////////////////////////////////////////////////////////////////
ARM_PricingModel* ARM_CRALocalCalculator::GetSFRMModel(void) const
{
	string ccy = string(itsCcy.GetCcyName());

	//GET MULTI ASSET MODEL
	ARM_MultiAssetsModel*  multiAssetsModel = NULL;
	multiAssetsModel = dynamic_cast<ARM_MultiAssetsModel*>(&*(GetPricingModel()));

	//SFRM MODEL
	ARM_PricingModel* sfrmModel = dynamic_cast<ARM_PricingModel*> (&*(*multiAssetsModel->GetModelMap())[ccy]->Model());

	return sfrmModel;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRALocalCalculator
///	Routine: Calibrate
///	Returns: void
///	Action : cfalibrate the model 
/////////////////////////////////////////////////////////////////
void ARM_CRALocalCalculator::Calibrate()
{
	try
	{
		//SFRM Calibration
		ARM_CRACalculator::Calibrate();

		//Local Model Calibration
		CalibrateLocalModel();
	}
	catch (Exception& x)
	{
	    x.DebugPrint();

		throw x;
	}
	catch (...)
	{
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in CRALocalCalculator::Calibrate" );
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRALocalCalculator
///	Routine: Price
///	Returns: a double
///	Action : price the Bermuda Swaption deal.
/////////////////////////////////////////////////////////////////
double ARM_CRALocalCalculator::Price()
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

//ARM_Timer timer;
//timer.ClockStartTime();

		genPricer->Price();

//timer.ClockEndTime();
//FILE* f=fopen("c:\\temp\\dumpCRA2Timer.txt","a");
//fprintf(f,"Price Duration = %10.5lf ms\n",timer.GetDuration()*1000.0);
//fclose(f);

		int NbProductsToPrice;
		
		//Local CapFloor Model
		if (itsLocalModelType == LocalDownUp)
		{
			NbProductsToPrice = MIN(localCfNbProductsToPrice, itsProductsToPrice.size());
			for (size_t i(0); i < NbProductsToPrice; ++i)
			{
				if (itsProductsToPrice[i])
				{
					price	= genPricer->GetPricerInfo()->GetContents( LocalCapFloorCRAColNamesTable[ LocalCapFloorCRAProductToPriceColumns[i] ] ).GetData("Price").GetDouble();
						
					if (i == localCfBermudaPrice)
					{
						itsBermuda1Price = price;
					}
					else if (i == localCfCorridorPrice)
					{
						itsCorridorPrice = price;
					}
					else if (i == localCfFundingPrice)
					{
						itsFundingPrice = price;
					}
					else if (i == localCfFwdCorridorPrice)
					{
						itsFwdCorridorPrice = price;
					}
					else if (i == localCfFwdFundingPrice)
					{
						itsFwdFundingPrice = price;
					}
					else if (i == localCfCorridorLegPrices)
					{
						itsCorridorLegPrices = std::vector<double>( * (genPricer->GetPricerInfo()->GetContents( LocalCapFloorCRAColNamesTable[ LocalCapFloorCRAProductToPriceColumns[i] ] ).GetData("IntermediatePrices").GetVector()) );
					}
				}
			}
			if(itsIsExerciseProbas)
			{
				double numValue0 = GetPricingModel()->GetNumeraire()->GetValue0();
				size_t nbExers = genSec->GetDealDescription().GetRowsNb()-1;
				itsExerciseProbas.resize(nbExers);
				size_t nbProbas = CC_Min(LocalCapFloorCRAExerciseProbaNb,nbExers-itsExerciseProbaOffset);
				for(i=0;i<itsExerciseProbaOffset;++i)
					itsExerciseProbas[i]=0.0;
				double prevProba = 0.0,proba;
				for(;i<itsExerciseProbaOffset+nbProbas;++i)
				{
					proba = genPricer->GetPricerInfo()->GetContents( LocalCapFloorCRAExerciseProbaNamesTable[i-itsExerciseProbaOffset] ).GetData("Price").GetDouble()
							/ numValue0;
					itsExerciseProbas[i] = proba - prevProba;
					prevProba = proba;
				}
				for(;i<nbExers;++i)
					itsExerciseProbas[i]=0.0;
			}
		}
		//Local Floor Model
		else
		{
			NbProductsToPrice = MIN(localFNbProductsToPrice, itsProductsToPrice.size());
			for (size_t i(0); i < NbProductsToPrice; ++i)
			{
				if (itsProductsToPrice[i])
				{
					price = genPricer->GetPricerInfo()->GetContents( LocalFloorCRAColNamesTable[ LocalFloorCRAProductToPriceColumns[i] ] ).GetData("Price").GetDouble();
						
					if (i == localFBermudaPrice)
					{
						itsBermuda1Price = price;
					}
					else if (i == localFCorridorPrice)
					{
						itsCorridorPrice = price;
					}
					else if (i == localFFundingPrice)
					{
						itsFundingPrice = price;
					}
					else if (i == localFFwdCorridorPrice)
					{
						itsFwdCorridorPrice = price;
					}
					else if (i == localFFwdFundingPrice)
					{
						itsFwdFundingPrice = price;
					}
				}
			}
		}

		itsHasBeenPriced = true;

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
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in CraLocalCalculator::Price" );
	}
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRALocalCalculator
///	Routine: ComputePricingData
///	Returns: a ARM_MultiTypeDict
///	Action : get pricing data of the calculaltor
/////////////////////////////////////////////////////////////////
void ARM_CRALocalCalculator::ComputePricingData() const
{
	if (!itsHasBeenPriced)
		const_cast<ARM_CRALocalCalculator*>(this)->PriceAndTimeIt();

	GetPricingData()[ "Bermuda" ]		= itsBermuda1Price;
	GetPricingData()[ "Corridor" ]		= itsCorridorPrice;
	GetPricingData()[ "Funding" ]		= itsFundingPrice;	
	GetPricingData()[ "FwdCorridor" ]	= itsFwdCorridorPrice;
	GetPricingData()[ "FwdFunding" ]	= itsFwdFundingPrice;	

	/// Prices of each residual corridor leg
	GetPricingData()[ "CorridorLeg" ]		= ARM_VectorPtr(new std::vector<double>(itsCorridorLegPrices));
	
	/// Exercice probabilities
	if (itsIsExerciseProbas)
		GetPricingData()[ "ExerciseProbas" ] = ARM_VectorPtr(new std::vector<double>(itsExerciseProbas));

}

////////////////////////////////////////////////////
///	Class   : ARM_CRALocalCalculator
///	Routines: toString
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
string ARM_CRALocalCalculator::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream craData;

    return craData.str();
}

////////////////////////////////////////////////////
///	Class   : ARM_CRALocalCalculator
///	Routines: GeneralDataToString
///	Returns :
///	Action  : Construct information string for
///           general product data
////////////////////////////////////////////////////
string ARM_CRALocalCalculator::GeneralDataToString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream generalCraData;

	//General Datas
	generalCraData << "GENERAL DATAS"	<<  endl;
	generalCraData << "Currency\t\t: "	<<	itsCcy.GetCcyName() << endl;
	generalCraData << "StartDate\t\t: " <<  itsStartDate.toString() << endl;
    generalCraData << "EndDate\t\t: "	<<	itsEndDate.toString() << endl;
    generalCraData << "Pay/Rec\t\t: "	<<  ARM_ParamView::GetMappingName(S_RECEIVE_PAY, itsPayRec) << endl;
	generalCraData << "\n";
	
	//Call Datas
	generalCraData << "CALL DATAS"		<< endl;
	generalCraData << "\nCall Frequency\t\t: "  << ARM_ArgConvReverse_StdFrequency.GetString(itsCallFreq) << endl;
   	generalCraData << "Call Notice\t\t:  "		<< itsCallNotice << "\n";
	generalCraData << "Call Calendar\t\t: "		<< itsCallCal << "\n";
	generalCraData << "\n";

	//Fund Datas
	generalCraData << "FUND DATAS"		<< endl;
    generalCraData << "Fund Frequency\t\t: " << ARM_ArgConvReverse_StdFrequency.GetString(itsFundFreq) << endl;
	generalCraData << "Fund Day Count\t\t: " << ARM_ArgConvReverse_DayCount.GetString(itsFundDayCount)  << endl;
	generalCraData << "\n";

	//Coupon Datas
	generalCraData << "EXOTIC COUPON DATAS" << endl;
    generalCraData << "Cpn Pay Freq\t\t: "			<< ARM_ArgConvReverse_StdFrequency.GetString(itsCpnPayFreq) << endl;
   	generalCraData << "Cpn Reset Calendar\t\t: "	<< itsCpnResetCal << "\n";
	generalCraData << "Cpn Pay Calendar\t\t: "		<< itsCpnPayCal << "\n";

	generalCraData << "Boosted Index Type\t\t: "	<< ARM_ArgConvReverse_IndexClass.GetString(itsBoostedIndexType) << endl;
	
	if (itsBoostedIndexType != K_FIXED)
	{
		generalCraData << "Boosted Var Term\t\t: "		<< itsBoostedVarTerm << endl;
		generalCraData << "Boosted Reset Gap\t\t: "		<< itsBoostedResetGap << endl;
		generalCraData << "Boosted Reset Timing\t\t: "	<< ARM_ArgConvReverse_Timing.GetString(itsBoostedResetTiming) << endl;
		generalCraData << "Boosted Day Count\t\t: "		<< ARM_ArgConvReverse_DayCount.GetString(itsBoostedDayCount)  << endl;
		generalCraData << "Boosted Adj Rule\t\t: "		<< ARM_ArgConvReverse_FwdRules.GetString(itsBoostedAdjRule)  << endl;
		generalCraData << "Boosted Rule\t\t: "			<< ARM_ArgConvReverse_InterestRules.GetString(itsBoostedIntRule )  << endl;
	}

	generalCraData << "Cpn Reset Freq\t\t: "			<< ARM_ArgConvReverse_StdFrequency.GetString(itsCpnResetFreq) << endl;
    generalCraData << "Cpn Reset Timing\t\t: "			<< ARM_ArgConvReverse_Timing.GetString(itsCpnResetTiming) << endl;
	generalCraData << "Cpn Reset Gap\t\t: "			<< itsCpnResetGap << "\n";

	generalCraData << "Ref Index Type 1\t\t: "			<< ARM_ArgConvReverse_IndexClass.GetString(itsRefIndexType1) << endl;
	generalCraData << "Ref Term 1\t\t: "				<< itsRefTerm1 << "\n";
	generalCraData << "Ref Day Count 1\t\t: "			<< ARM_ArgConvReverse_DayCount.GetString(itsRefDayCount1)  << endl;
	generalCraData << "Ref Coeff 1\t\t: "				<< itsRefCoeff1 << "\n";


	generalCraData << endl;

    return generalCraData.str();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRALocalCalculator
///	Routine: View
///	Returns: 
///	Action : .
/////////////////////////////////////////////////////////////////
void ARM_CRALocalCalculator::View(char* id, FILE* ficOut) const
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
    fprintf(fOut,"\n\n ==================================================================== \n");
	fprintf(fOut,"\n\n =======>  CALLABLE RANGE ACCRUAL LOCAL CALCULATOR <================= \n");
    fprintf(fOut,"\n\n ==================================================================== \n");

	fprintf(fOut,"%s",GeneralDataToString().c_str());
	fprintf(fOut,"%s",DealDesDataToString().c_str());

	CC_Ostringstream dataLocal;
	if (itsLocalModelType == LocalDown)
	{
		fprintf(fOut,"\n\n ============>    LOCAL FLOOR CALIBRATION    <=========== \n");
	}
	else if (itsLocalModelType == LocalUp)
	{
		fprintf(fOut,"\n\n ============>    LOCAL CAP CALIBRATION    <=========== \n");
	}
	else
	{
		fprintf(fOut,"\n\n ============>    LOCAL CAP FLOOR CALIBRATION    <=========== \n");	
	}

	CC_Ostringstream localPf;
	fprintf(fOut,"\n\n ============>    LOCAL CORRIDOR PORTFOLIO     <=========== \n");
	localPf << itsLocalPortfolio->toString() << "\n";
	fprintf(fOut,"%s",localPf.str().c_str());
	fprintf(fOut,"\n\n ============>  EOF LOCAL CORRIDOR PORTFOLIO   <=========== \n");

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

	if ( ficOut == NULL )
       fclose(fOut);
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
