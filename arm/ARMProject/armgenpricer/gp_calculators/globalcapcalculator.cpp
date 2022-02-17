
/*!
 /** Copyright (c) CDC IXIS CM July 2003 Paris
 *	\file GlobalCapcalculator.cpp
 *  \brief file for the Global Cap Calculator
 *	\author  P. LAM & H.ESSEFI
 *	\version 1.0
 *	\date Feb 2006
 */

/// this header comes first as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpcalculators/globalcapcalculator.h"
#include "gpcalculators/argconvdefault.h"

/// gpbase
#include "gpbase/argconvdefault.h"
#include "gpbase/datestripcombiner.h"
#include "gpbase/autocleaner.h"
#include "gpbase/singleton.h"
#include "gpbase/datestrip.h"
#include "gpbase/curveconvert.h"
#include "gpbase/curve.h"
#include "gpbase/gplinalgconvert.h"
#include "gpbase/stringmanip.h"

/// gpinfra
#include "gpinfra/argconvdefault.h"
#include "gpinfra/gensecurity.h"
#include "gpinfra/dealdescription.h"
#include "gpinfra/mktdatamanagerrep.h"
#include "gpinfra/pricingadviser.h"
#include "gpinfra/genpricer.h"
#include "gpinfra/numeraire.h"
#include "gpinfra/numerairefactory.h"
#include "gpinfra/pricerinfo.h"
#include "gpinfra/cstmanager.h"
#include "gpinfra/curvemodelparam.h"

/// gpcalib
#include "gpcalib/calibmethod.h"
#include "gpcalib/vanillasecuritydensity.h"

/// gpmodels
#include "gpmodels/MarkovFunctional.h"
#include "gpmodels/ModelParamsMF.h"

/// gpnumlib
#include "gpnumlib/randomgenfactory.h"
#include "gpnumlib/antitheticgen.h"
#include "gpnumlib/argconvdefault.h"

/// gpnummethods
#include "gpnummethods/scheduler.h"
#include "gpnummethods/sampler.h"
#include "gpnummethods/truncator.h"
#include "gpnummethods/schedulerfactory.h"
#include "gpnummethods/samplerfactory.h"
#include "gpnummethods/pdemethod.h"
#include "gpnummethods/argconvdefault.h"
#include "gpnummethods/amcmethod.h"

/// kernel
#include "inst/swapleg.h"
#include "inst/fixleg.h"
#include "inst/portfolio.h"
#include "inst/barrier.h"
#include "inst/globalcap.h"

CC_BEGIN_NAMESPACE( ARM )

const int NB_FACTORS				= 1;

/// Default MDM key names
const string YC_KEY_NAME            = "YC_";
const string CFMODEL_KEY_NAME       = "CFMOD_";
const string MRS_KEY_NAME           = "MRS_";

/// Cap schedule 
const unsigned int GLOBALCAP_SCHED	= 0;

const int DEAL_NBCOLUMNS			= 21;

const string ARM_GlobalCapCalculator::GlobalCapColNamesTable [] =
{  
	"ResetDate",
	"StartDate",
	"EndDate",
	"PayDate",
	"Notional",
	"Spread",
	"Fixed",
	"CapLev",
	"FundLev", 
	"Strike",
	"Basis",
	"Index",
	"DF",
	"FundRate",
	"CouponRate",
	"Funding",
	"SumCoup",
	"Coupon",
	"Swap",
	"GlobalCap",
	"Product",
};


const int ARM_GlobalCapCalculator::GlobalCapProductToPriceColumns [] =
{
	Product,
	Funding,
	Coupon,
	Swap,
	GlobalCap,
	Index,
	FundRate,
	SumCoup,
	DF,
};


/*****************************************************************************
/ Class  : ARM_GlobalCapCalculator
/ Routine: Constructor 
/ Returns: void 
/ Action : builds the object 
******************************************************************************/
ARM_GlobalCapCalculator::ARM_GlobalCapCalculator( ARM_Currency& ccy,
												  ARM_Date&	startDate,
												  ARM_Date&	endDate,
												  int payRec,
												  ARM_ReferenceValue& notional,
												  ARM_INDEX_TYPE indexType,
												  int dayCount,
												  int payFreq,
												  int resetGap,
												  int payGap,
												  int resetTiming,
												  int payTiming,
												  string resetCal,
												  string payCal,
												  int adjRule,
												  int intRule,
												  ARM_ReferenceValue& fundLev,
												  ARM_ReferenceValue& capLev,
												  ARM_ReferenceValue& capFixed,
												  ARM_ReferenceValue& capStrike,
												  ARM_ReferenceValue& capSpread,
												  ARM_Vector* globalCapParams,
												  ARM_ReferenceValue* pastFixings,
												  int nbSteps,
												  vector<string> randGenerator,
												  int samplerType,
												  ARM_Vector* calibParams,
												  ARM_StringVector& mdmKeys,
												  const ARM_MarketData_ManagerRep& mktDataManager,
												  const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice)
: ARM_GenCalculator(mktDataManager),
itsCcy(ccy),
itsStartDate(startDate),
itsEndDate(endDate),
itsPayRec(payRec),
itsNotional(notional),
itsIndexType(indexType),
itsDayCount(dayCount),
itsPayFreq(payFreq),
itsResetGap(resetGap),
itsPayGap(payGap),
itsResetTiming(resetTiming),
itsPayTiming(payTiming),
itsResetCal(resetCal),
itsPayCal(payCal),
itsAdjRule(adjRule),
itsIntRule(intRule),
itsFundLev(fundLev),
itsCapLev(capLev),
itsCapFixed(capFixed),
itsCapStrike(capStrike),
itsCapSpread(capSpread),
itsGlobalCapParams(globalCapParams),
itsPastFixings(pastFixings),
itsNbSteps(nbSteps),
itsRandGenerator(randGenerator),
itsSamplerType(samplerType),
itsCalibParams(calibParams),
itsMCMethod(NULL),
itsProductsToPrice(productsToPrice),
itsHasBeenPriced(false),
itsFundingPrice(0.0),
itsGlobalCapPrice(0.0),
itsCouponPrice(0.0),
itsSwapPrice(0.0),
itsProductPrice(0.0),
itsSumCoupPrice(0.0),
itsDFPrice(0.0),
itsFundRatePrice(0.0),
itsIndexPrice(0.0)
{
	SetCurrencyUnit(&ccy); 

	SetKeys(mdmKeys);

	CheckData();

	ARM_StringVector pricedColumns = PricedColumnNames();

	ARM_CstManagerPtr cstManagerPtr = CreateCstManager();

    CreateAndSetDealDescriptionAndTimeIt("", pricedColumns, cstManagerPtr);

	CreateAndSetModel();

	CreateAndSetCalibration(); 
}

/*****************************************************************************
/ Class  : ARM_GlobalCapCalculator 
/ Routine: Constructor from Cap
/ Returns: void 
/ Action : builds the object 
******************************************************************************/
ARM_GlobalCapCalculator::ARM_GlobalCapCalculator( ARM_GlobalCap* globalCap,
												  ARM_ReferenceValue& fundLev,
												  ARM_ReferenceValue& capLev,
												  int nbSteps,
												  vector<string> randGenerator,
												  int samplerType,
												  ARM_Vector* calibParams,
												  ARM_StringVector& mdmKeys,
												  const ARM_MarketData_ManagerRep& mktDataManager,
												  const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice)
: ARM_GenCalculator(mktDataManager),
itsFundLev(fundLev),
itsCapLev(capLev),
itsNbSteps(nbSteps),
itsRandGenerator(randGenerator),
itsSamplerType(samplerType),
itsCalibParams(calibParams),
itsProductsToPrice(productsToPrice),
itsMCMethod(NULL),
itsHasBeenPriced(false),
itsFundingPrice(0.0),
itsGlobalCapPrice(0.0),
itsCouponPrice(0.0),
itsSwapPrice(0.0),
itsProductPrice(0.0),
itsSumCoupPrice(0.0),
itsDFPrice(0.0),
itsFundRatePrice(0.0),
itsIndexPrice(0.0),
itsPastFixings(NULL)
{
	SetName(ARM_GLOBALCAP_CALCULATOR);

	ARM_SwapLeg* fundingLeg = globalCap->GetSwapLeg();
	ARM_Currency* ccy = globalCap->GetCurrencyUnit();
	SetCcy(*ccy);
	SetStartDate(fundingLeg->GetStartDateNA());
	SetEndDate(fundingLeg->GetEndDateNA());
	SetPayRec(globalCap->GetPorS());
	SetNotional(*(globalCap->GetAmount()));

	int idxtyp = ccy->GetVanillaIndexType();
	SetIndexType(fundingLeg->GetIRIndex()->GetIndexType());
	SetDayCount(fundingLeg->GetDayCount()); 
	SetPayFreq(fundingLeg->GetPaymentFreq());
	SetResetGap(fundingLeg->GetIRIndex()->GetResetGap());
	SetPayGap(fundingLeg->GetIRIndex()->GetPayGap());
	SetResetTiming(fundingLeg->GetIRIndex()->GetResetTiming());
	SetPayTiming(fundingLeg->GetIRIndex()->GetPayTiming());
	SetResetCal(fundingLeg->GetResetCalName());
	SetPayCal(fundingLeg->GetPayCalName());
	SetAdjRule(fundingLeg->GetDecompPricingFlag());
	SetIntRule(fundingLeg->GetIRIndex()->GetIntRule());
	SetCapFixed(*(globalCap->GetFixedRates()));
	SetCapStrike(*(globalCap->GetBarriers()));
	SetCapSpread(*(globalCap->GetSpreads()));

	// Global Cap Params
	ARM_Vector* globalCapParams = new ARM_Vector(3);
	
	globalCapParams->InitElt(0, (double) globalCap->GetSwapLeg()->GetResetDates()->size());
	globalCapParams->InitElt(1, globalCap->GetRatio());
	globalCapParams->InitElt(2, globalCap->GetStrike() /100);

	SetGlobalCapParams(globalCapParams);

	delete globalCapParams;
	globalCapParams = NULL;

	SetKeys(mdmKeys);

	CheckData();

	ARM_StringVector pricedColumns = PricedColumnNames();

	ARM_CstManagerPtr cstManagerPtr = CreateCstManager();

	CreateAndSetDealDescriptionAndTimeIt("", pricedColumns, cstManagerPtr);

	CreateAndSetModel();

	CreateAndSetCalibration();
}


/*****************************************************************************
/ Class  : ARM_GlobalCapCalculator 
/ Routine: Constructor from Cap (requires Init after) 
/ Returns: void 
/ Action : builds the object 
******************************************************************************/

ARM_GlobalCapCalculator::ARM_GlobalCapCalculator( ARM_Date& asOfDate,
												  ARM_GlobalCap* globalCap,
												  ARM_ReferenceValue& fundLev,
												  ARM_ReferenceValue& capLev)
:ARM_GenCalculator(asOfDate),
itsFundLev(fundLev),
itsCapLev(capLev),
itsNbSteps(0),
itsRandGenerator(NULL),
itsSamplerType(0),
itsCalibParams(NULL),
itsProductsToPrice(NULL),
itsMCMethod(NULL),
itsHasBeenPriced(false),
itsFundingPrice(0.0),
itsGlobalCapPrice(0.0),
itsCouponPrice(0.0),
itsSwapPrice(0.0),
itsProductPrice(0.0),
itsSumCoupPrice(0.0),
itsDFPrice(0.0),
itsFundRatePrice(0.0),
itsIndexPrice(0.0),
itsPastFixings(NULL)
{
	SetName(ARM_GLOBALCAP_CALCULATOR);

	ARM_SwapLeg* fundingLeg = globalCap->GetSwapLeg();
	ARM_Currency* ccy = globalCap->GetCurrencyUnit();
	SetCcy(*ccy);

	ARM_Vector* myResetDates = fundingLeg->GetResetDates();
	double startDate;
	bool fundStartDate (false);

	for (int i = 0; i < myResetDates->GetSize(); i++)
	{
		if (myResetDates->Elt(i) > asOfDate.GetJulian())
		{
			startDate = fundingLeg->GetFlowStartDates()->Elt(i);
			fundStartDate = true;
			break;
		}
	}

	if (!fundStartDate)
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in GlobalCapCalculator: occured in the past" );

	SetStartDate(ARM_Date(startDate));
	SetEndDate(fundingLeg->GetEndDateNA());
	SetPayRec(globalCap->GetPorS());
	SetNotional(*(globalCap->GetAmount()));

	int idxtyp = ccy->GetVanillaIndexType();
	SetIndexType(fundingLeg->GetIRIndex()->GetIndexType());
	SetDayCount(fundingLeg->GetDayCount()); 
	SetPayFreq(fundingLeg->GetPaymentFreq());
	SetResetGap(fundingLeg->GetIRIndex()->GetResetGap());
	SetPayGap(fundingLeg->GetIRIndex()->GetPayGap());
	SetResetTiming(fundingLeg->GetIRIndex()->GetResetTiming());
	SetPayTiming(fundingLeg->GetIRIndex()->GetPayTiming());
	SetResetCal(fundingLeg->GetResetCalName());
	SetPayCal(fundingLeg->GetPayCalName());
	SetAdjRule(fundingLeg->GetDecompPricingFlag());
	SetIntRule(fundingLeg->GetIRIndex()->GetIntRule());
	SetCapFixed(*(globalCap->GetFixedRates()));
	SetCapStrike(*(globalCap->GetBarriers()));
	SetCapSpread(*(globalCap->GetSpreads()));
	SetPastFixings(globalCap->GetPastFixings());

	ARM_Vector* globalCapParams = new ARM_Vector(3);

	globalCapParams->InitElt(0, (double) globalCap->GetSwapLeg()->GetResetDates()->size());
	globalCapParams->InitElt(1, globalCap->GetRatio());
	globalCapParams->InitElt(2, globalCap->GetStrike() /100);

	SetGlobalCapParams(globalCapParams);
	delete globalCapParams;
	globalCapParams = NULL;

	double lastPastFixingDate;
	double nextFixingDate;

	if (itsPastFixings)
	{
		lastPastFixingDate = myResetDates->Elt(globalCap->GetPastFixings()->GetDiscreteDates()->GetSize() - 1);
		nextFixingDate = myResetDates->Elt(globalCap->GetPastFixings()->GetDiscreteDates()->GetSize());
	}
	else
	{
		lastPastFixingDate = ARM_Date(ARM_DEFAULT_DATE).GetJulian();
		nextFixingDate = myResetDates->Elt(0);
	}

	if ( lastPastFixingDate < asOfDate.GetJulian() && nextFixingDate < asOfDate.GetJulian() )	// Il manque des fixings
	{
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in GlobalCapCalculator: past fixings missing" );
	}
	else if ( nextFixingDate == asOfDate.GetJulian() )
	{
		// Le fixing a lieu en date asOf mais n'a pas encore été renseigné		
		ARM_Vector* pastFixingDates = NULL;
		ARM_Vector* pastFixingValues = NULL;

		if (itsPastFixings)
		{
			pastFixingDates = itsPastFixings->GetDiscreteDates();
			pastFixingValues = itsPastFixings->GetDiscreteValues();
			
			pastFixingDates->push_back(asOfDate.GetJulian());
			pastFixingValues->push_back(-1);
		}
		else
		{
			pastFixingDates = new ARM_Vector();
			pastFixingValues = new ARM_Vector();

			pastFixingDates->push_back(asOfDate.GetJulian());
			pastFixingValues->push_back(-1);

			ARM_ReferenceValue* tmpFixings = new ARM_ReferenceValue(pastFixingDates, pastFixingValues);
			SetPastFixings(tmpFixings);
			
			if (tmpFixings)
				delete tmpFixings;
			tmpFixings = NULL;
		}
	}
}

/*****************************************************************************
/ Class  : InitGlobalCapFromSummit
/ Routine: 
/ Returns: void 
/ Action : Get Data From Summit 
******************************************************************************/

void ARM_GlobalCapCalculator::InitGlobalCapFromSummit(ARM_ZeroCurve* zc, 
													  ARM_VolCurve* capVol,
													  ARM_VolCurve* rhoCap,
													  ARM_VolCurve* nuCap,
													  ARM_VolCurve* betaCap,
													  double mrs,
													  ARM_Vector* calibParams,
													  int nbSteps,
													  vector<string> randGenerator,
													  int samplerType,
													  const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice)
{
	SetCalibParams(calibParams);
	itsNbSteps = nbSteps;
	itsRandGenerator = randGenerator;
	itsSamplerType = samplerType;
	itsProductsToPrice = productsToPrice;

	// Définit les données à importer de Summit
	
	ARM_StringVector keys(NbKeys);
	string ccyName(itsCcy.GetCcyName());

	keys[YcKey]			= YC_KEY_NAME + ccyName;
    keys[CfModelKey]	= CFMODEL_KEY_NAME	+ ccyName;
    keys[MrsKey]		= MRS_KEY_NAME	+ ccyName;

	SetKeys(keys);

	ARM_MarketData_ManagerRep *myMktDataManager = &*GetMktDataManager();	

	ARM_Date asof = myMktDataManager->GetAsOfDate();
	
	// Initialise le MktDataManager
	//InitializeMktDataManagerOnly(*myMktDataManager);
	
	if (itsPastFixings)
	{
		int length = itsPastFixings->GetDiscreteDates()->GetSize();

		if ( itsPastFixings->GetDiscreteDates()->Elt(length - 1) == asof.GetJulian() && itsPastFixings->GetDiscreteValues()->Elt(length - 1) == -1)
		{
			itsPastFixings->GetDiscreteValues()->Elt(length - 1) = zc->DiscountYield(FromIndexTypeToTerm(itsIndexType))/100;
		}
	}

	// Cap Model
	ARM_BSSmiledModel* capSabrModel = NULL;
	ARM_VolCurve* ATMCapVol = (capVol->GetName() == ARM_VOL_CUBE) ? ((ARM_VolCube*)capVol)->GetATMVol() : capVol;
	
	int methodType;
	if (betaCap)
	{
		if (betaCap->IsEqualToOne())
			methodType = K_SABR_ARITH;
		else
			methodType = K_SABR_IMPLNVOL;
	}
	
	capSabrModel = new ARM_BSSmiledModel(	asof, 
											0.0, 
											zc, 
											zc, 
											ATMCapVol, 
											K_YIELD, 
											rhoCap, 
											nuCap, 
											methodType, 
											betaCap);
	
	std::vector<double> breakPointTimes(1, 0.0);
	std::vector<double> values(1, mrs);
	ARM_CurveModelParam* mrsCurve = new ARM_CurveModelParam(ARM_ModelParamType::MeanReversion, &values, &breakPointTimes);

	// Importe les données dans le Market Data Manager
	myMktDataManager->RegisterData(keys[YcKey], zc);
	myMktDataManager->RegisterData(keys[CfModelKey], capSabrModel);
	myMktDataManager->RegisterData(keys[MrsKey], mrsCurve);

	delete capSabrModel;
	capSabrModel = NULL;

	delete mrsCurve;
	mrsCurve = NULL;

	CheckData();

	ARM_StringVector pricedColumns = PricedColumnNames();
	ARM_CstManagerPtr cstManagerPtr = CreateCstManager();
	CreateAndSetDealDescriptionAndTimeIt("", pricedColumns, cstManagerPtr);
	CreateAndSetModel();
	CreateAndSetCalibration();
}


////////////////////////////////////////////////////
///	Class  : ARM_GlobalCapCalculator
///	Routine: CopyNoCleanUp
///	Returns: void
///	Action : Arguments copy
////////////////////////////////////////////////////
void ARM_GlobalCapCalculator::CopyNoCleanUp(const ARM_GlobalCapCalculator& rhs)
{
	itsCcy						= rhs.itsCcy;
	itsStartDate				= rhs.itsStartDate;
	itsEndDate					= rhs.itsEndDate;
	itsPayRec					= rhs.itsPayRec;
	itsNotional					= rhs.itsNotional;
	itsIndexType				= rhs.itsIndexType;
	itsDayCount					= rhs.itsDayCount;
	itsPayFreq					= rhs.itsPayFreq;
	itsResetGap					= rhs.itsResetGap;
	itsPayGap					= rhs.itsPayGap;
	itsResetTiming				= rhs.itsResetTiming;
	itsPayTiming				= rhs.itsPayTiming;
	itsResetCal					= rhs.itsResetCal;
	itsPayCal					= rhs.itsPayCal;
	itsAdjRule					= rhs.itsAdjRule;
	itsIntRule					= rhs.itsIntRule;
	itsFundLev					= rhs.itsFundLev;
	itsCapLev					= rhs.itsCapLev;
	itsCapFixed					= rhs.itsCapFixed;
	itsCapStrike				= rhs.itsCapStrike;
	itsCapSpread				= rhs.itsCapSpread;
	itsGlobalCapParams			= (rhs.itsGlobalCapParams ? (ARM_Vector *) rhs.itsGlobalCapParams->Clone() : NULL);
	itsNbSteps					= rhs.itsNbSteps;
	itsRandGenerator			= rhs.itsRandGenerator;
	itsSamplerType				= rhs.itsSamplerType;
	itsCalibParams				= (rhs.itsCalibParams ? (ARM_Vector *) rhs.itsCalibParams->Clone() : NULL);
	itsMCMethod					= rhs.itsMCMethod;
	itsProductsToPrice			= rhs.itsProductsToPrice;
	itsHasBeenPriced			= rhs.itsHasBeenPriced;
	itsFundingPrice				= rhs.itsFundingPrice;
	itsGlobalCapPrice			= rhs.itsGlobalCapPrice;
	itsCouponPrice				= rhs.itsCouponPrice;
	itsSwapPrice				= rhs.itsSwapPrice;
	itsProductPrice				= rhs.itsProductPrice;
	itsSumCoupPrice				= rhs.itsSumCoupPrice;
	itsDFPrice					= rhs.itsDFPrice;
	itsFundRatePrice			= rhs.itsFundRatePrice;
	itsIndexPrice				= rhs.itsIndexPrice;
	itsExerDateStrip			= rhs.itsExerDateStrip;

	if (rhs.itsPastFixings)
	{
		ARM_Vector* pastDates = (ARM_Vector*) rhs.itsPastFixings->GetDiscreteDates()->Clone();
		ARM_Vector* pastValues = (ARM_Vector*) rhs.itsPastFixings->GetDiscreteValues()->Clone();

		ARM_ReferenceValue pastFixings(pastDates, pastValues);
		SetPastFixings(&pastFixings);
	}
	else
	{
		itsPastFixings = NULL;
	}
}

/*****************************************************************************
/ Class  : ARM_GlobalCapCalculator
/ Routines: Clone
/ Returns :
/ Action  : Standard ARM object support
/           Call copy constructor
******************************************************************************/
ARM_Object* ARM_GlobalCapCalculator::Clone() const
{
	return(new ARM_GlobalCapCalculator(*this));
}

/*****************************************************************************
/ Class  : ARM_GlobalCapCalculator
/ Routines: 
/ Returns :
/ Action  : copy constructor
******************************************************************************/
ARM_GlobalCapCalculator::ARM_GlobalCapCalculator(const ARM_GlobalCapCalculator& rhs)
:ARM_GenCalculator(rhs)
{
	CopyNoCleanUp(rhs);
}

/*****************************************************************************
/ Class  : ARM_GlobalCapCalculator
/ Routines: operator=
/ Returns :
/ Action  : Call copy constructor
******************************************************************************/
ARM_GlobalCapCalculator& ARM_GlobalCapCalculator::operator=(const ARM_GlobalCapCalculator& rhs)
{
	if( this != & rhs )
	{
		ARM_GenCalculator::operator=(rhs);
        CleanUp();
        CopyNoCleanUp(rhs);
	}
	return *this;
}

/*****************************************************************************
/ Class  : ARM_GlobalCapCalculator
/ Routine: Destructor
/ Returns: void
/ Action : destroys the object
******************************************************************************/
ARM_GlobalCapCalculator::~ARM_GlobalCapCalculator()
{
	CleanUp();
}

/*****************************************************************************
/ Class  : ARM_GlobalCapCalculator
/ Routine: CleanUp
/ Returns: void
/ Action : delete pointors
******************************************************************************/
void ARM_GlobalCapCalculator::CleanUp()
{
	if (itsCalibParams)
	{
		delete itsCalibParams;
		itsCalibParams = NULL;
	}
	if (itsGlobalCapParams)
	{
		delete itsGlobalCapParams;
		itsGlobalCapParams = NULL;
	}
	if (itsPastFixings)
	{
		delete itsPastFixings;
		itsPastFixings = NULL;
	}
}

/*****************************************************************************
/ Class  : ARM_GlobalCapCalculator
/ Routine: CreateCstManager
/ Returns: void
/ Action : create the const manager (static data).
******************************************************************************/
ARM_CstManagerPtr ARM_GlobalCapCalculator::CreateCstManager()
{
	try
	{
		double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

		vector<string>	cstNames; 

		cstNames.push_back("NotionalCst");
		cstNames.push_back("FundLevCst");
		cstNames.push_back("CapLevCst");
		cstNames.push_back("CapFixedCst");
		cstNames.push_back("CapStrikeCst");
		cstNames.push_back("CapSpreadCst");

		vector<ARM_GramFctorArg> cstVector;
		
		//Notional
		ARM_Curve* tmpNotional = RefValueToCurve(itsNotional, asOfDate);
		cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpNotional))));

		//FundLev
		ARM_Curve* tmpFundLev = RefValueToCurve(itsFundLev, asOfDate);
		cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpFundLev))));

		//CapLev
		ARM_Curve* tmpCapLev = RefValueToCurve(itsCapLev, asOfDate);
		cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpCapLev))));

		//CapFixed
		ARM_Curve* tmpCapFixed = RefValueToCurve(itsCapFixed, asOfDate);
		cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpCapFixed))));

		//CapStrike
		ARM_Curve* tmpCapStrike = RefValueToCurve(itsCapStrike, asOfDate);
		cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpCapStrike))));

		//CapSpread
		ARM_Curve* tmpCapSpread = RefValueToCurve(itsCapSpread, asOfDate);
		cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpCapSpread))));

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
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in GlobalCapCalculator::CreateCstManager" );
	}
}


/*****************************************************************************
/ Class  : ARM_GlobalCapCalculator
/ Routine: CheckMktData
/ Returns: 
/ Action : Checks if market data are of good type
******************************************************************************/
void ARM_GlobalCapCalculator::CheckMktData()
{
    /// MdM datas checking
	ARM_ZeroCurve* ycCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));
    if(!ycCurve)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Yield Curve for key=" + GetKeys()[YcKey] + " is expected in the Market Data Manager");

	ARM_BSModel* cfBSModel = dynamic_cast< ARM_BSModel* >(GetMktDataManager()->GetData(GetKeys()[CfModelKey]) );
    if(!cfBSModel)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : cap/floor B&S smiled model for key=" + GetKeys()[CfModelKey] + " is expected in the Market Data Manager");

	ARM_ModelParam* mrsParam = dynamic_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[MrsKey]));
    if(!mrsParam || mrsParam->GetType() != ARM_ModelParamType::MeanReversion)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : MRS Param for key=" + GetKeys()[MrsKey] + " is expected in the Market Data Manager");
}

/*****************************************************************************
/ Class  : ARM_GlobalCapCalculator
/ Routine: CheckData
/ Returns: 
/ Action : For the moment, only checks market data
******************************************************************************/
void ARM_GlobalCapCalculator::CheckData()
{
	CheckMktData();
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_GlobalCapCalculator
///	Routine: CustomDatesStructure
///	Returns: ARM_DateStripCombiner
///	Action : create the list of all event dates of the Global Cap. 
///    	 The DateStripCombiner merges event dates of each leg
///		customized Schedule Not implemented yet.
/////////////////////////////////////////////////////////////////
ARM_DateStripCombiner ARM_GlobalCapCalculator::CustomDatesStructure(const ARM_DateStripVector& dateStrips) const
{
	return ARM_DateStripCombiner();
}

/*****************************************************************************
/ Class  : ARM_GlobalCapCalculator
/ Routine: DatesStructure
/ Returns: ARM_DateStripVector
/ Action : create the list of all event dates of the Global Cap. 
/    	 The DateStripCombiner merges event dates of each leg
******************************************************************************/
ARM_DateStripCombiner ARM_GlobalCapCalculator::DatesStructure() const
{
	try
	{
		ARM_DateStripVector SchedVect(1,NULL);
	
		double startDate = itsStartDate.GetJulian();
		double endDate = itsEndDate.GetJulian();

		ARM_DateStrip ProdSched(startDate, endDate, itsPayFreq, itsDayCount,
								itsResetCal.c_str(), itsAdjRule, itsIntRule, 
								K_SHORTSTART, itsResetGap, itsPayFreq, GETDEFAULTVALUE,
								itsPayCal.c_str(), itsResetTiming, itsPayTiming, true);

		SchedVect[GLOBALCAP_SCHED] = &ProdSched;

		itsExerDateStrip = ARM_DateStripPtr(static_cast<ARM_DateStrip*>(ProdSched.Clone()));
	
		return SchedVect;
	}
	catch (Exception& x)
	{
	    x.DebugPrint();

		throw x;
	}
	catch (...)
	{
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in GlobalCapCalculator::DatesStructure" );
	}
}

/*****************************************************************************
/ Class  : ARM_GlobalCapCalculator
/ Routine: MiddleRows
/ Returns: ARM_RowInfo
/ Action : create a row of a deal description
******************************************************************************/
ARM_RowInfo ARM_GlobalCapCalculator::MiddleRows(size_t eventIdx, 
                                          const ARM_DateStripCombiner& datesStructure) const
{
	try
	{
		//Useful data
		char* ccy	= GetCcy().GetCcyName();
		double endDate = GetEndDate().GetJulian();
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

		//Number of columns
		size_t descSize = DEAL_NBCOLUMNS;
	
		size_t nbEvent = datesStructure.GetDateStrip(GLOBALCAP_SCHED)->GetResetDates()->size();
		vector< string > rowDescVec(descSize);
		vector< ARM_GP_VALUE_TYPE > rowTypeVec(descSize, ARM_MISSING_TYPE); 
		
		InitPriceableColumns(rowDescVec,rowTypeVec);

		//RESET DATE
		double resetDate;
		resetDate = (*(datesStructure.GetDateStrip(GLOBALCAP_SCHED)->GetResetDates()))[eventIdx];

		CC_Ostringstream resetDateDesc;
		resetDateDesc << CC_NS(std,fixed) << resetDate;
		rowDescVec[ResetDate] = resetDateDesc.str();
		rowTypeVec[ResetDate] = ARM_DATE_TYPE;

		//START DATE
		double startDate =  (*(datesStructure.GetDateStrip(GLOBALCAP_SCHED)->GetFwdRateStartDates()))[eventIdx] ;
		CC_Ostringstream startDateDesc;
		startDateDesc << CC_NS(std,fixed) << startDate;
		rowDescVec[StartDate] = startDateDesc.str();
		rowTypeVec[StartDate] = ARM_DATE_TYPE;

		//PAY DATE
		double payDate =  (*(datesStructure.GetDateStrip(GLOBALCAP_SCHED)->GetPaymentDates()))[eventIdx] ;
		CC_Ostringstream payDateDesc;
		payDateDesc << CC_NS(std,fixed) << payDate;
		rowDescVec[PayDate] = payDateDesc.str();
		rowTypeVec[PayDate] = ARM_DATE_TYPE;
		
		//NEXT START DATE / NEXT CALL DATE
		double nextStartDate = (*(datesStructure.GetDateStrip(GLOBALCAP_SCHED)->GetFwdRateEndDates()))[eventIdx] ;
		CC_Ostringstream nextStartDateDesc;
		nextStartDateDesc << CC_NS(std,fixed) << nextStartDate;
		rowDescVec[EndDate] = nextStartDateDesc.str();
		rowTypeVec[EndDate] = ARM_DATE_TYPE;

		//NOTIONAL
		double notional = const_cast< ARM_ReferenceValue& >(itsNotional).Interpolate(startDate);
		CC_Ostringstream notionalDesc;
		notionalDesc << CC_NS(std,fixed) << notional;
		rowDescVec[Notional] = notionalDesc.str();
		rowTypeVec[Notional] = ARM_DOUBLE;

		//SPREAD
		double spread = const_cast< ARM_ReferenceValue& >(itsCapSpread).Interpolate(startDate);
		CC_Ostringstream spreadDesc;
		spreadDesc << CC_NS(std,fixed) << spread;
		rowDescVec[Spread] = spreadDesc.str();
		rowTypeVec[Spread] = ARM_DOUBLE;

		//FIXED
		double fixed = const_cast< ARM_ReferenceValue& >(itsCapFixed).Interpolate(startDate);
		CC_Ostringstream fixedDesc;
		fixedDesc << CC_NS(std,fixed) << fixed;
		rowDescVec[Fixed] = fixedDesc.str();
		rowTypeVec[Fixed] = ARM_DOUBLE;

		//CAPLEV
		double capLev = const_cast< ARM_ReferenceValue& >(itsCapLev).Interpolate(startDate);
		CC_Ostringstream capLevDesc;
		capLevDesc << CC_NS(std,fixed) << capLev;
		rowDescVec[CapLev] = capLevDesc.str();
		rowTypeVec[CapLev] = ARM_DOUBLE;

		//FUNDLEV
		double fundLev = const_cast< ARM_ReferenceValue& >(itsFundLev).Interpolate(startDate);
		CC_Ostringstream fundLevDesc;
		fundLevDesc << CC_NS(std,fixed) << fundLev;
		rowDescVec[FundLev] = fundLevDesc.str();
		rowTypeVec[FundLev] = ARM_DOUBLE;

		//STRIKE
		double strike = const_cast< ARM_ReferenceValue& >(itsCapStrike).Interpolate(startDate);
		CC_Ostringstream strikeDesc;
		strikeDesc << CC_NS(std,fixed) << strike;
		rowDescVec[Strike] = strikeDesc.str();
		rowTypeVec[Strike] = ARM_DOUBLE;

		//BASIS
		double basis = CountYears(itsDayCount, startDate, nextStartDate);
		CC_Ostringstream basisDesc;
		basisDesc << CC_NS(std,fixed) << basis;
		rowDescVec[Basis] = basisDesc.str();
		rowTypeVec[Basis] = ARM_DOUBLE;

		//INDEX
		CC_Ostringstream indexDesc;
		indexDesc << "LIBOR(" << ccy << ",";
		indexDesc << "StartDate[i],";
		indexDesc << YearTermToStringMatu(FromIndexTypeToTerm(itsIndexType)) << ")";
		rowDescVec[Index] = indexDesc.str();
		rowTypeVec[Index] = ARM_STRING;

		//DF
		CC_Ostringstream dfDesc;
		dfDesc << "DF(" << ccy << ",";
		dfDesc << "PayDate[i])";
		rowDescVec[DF] = dfDesc.str();
		rowTypeVec[DF] = ARM_STRING;

		//FUNDRATE
		CC_Ostringstream fundRateDesc;
		fundRateDesc << "Index[i]*Basis[i]*DF[i]";
		rowDescVec[FundRate] = fundRateDesc.str();
		rowTypeVec[FundRate] = ARM_STRING;

		//COUPONRATE
		CC_Ostringstream couponRateDesc;
		couponRateDesc << "CapLev[i]*Index[i]+if(Index[i]<Strike[i],Fixed[i],Index[i]+Spread[i])";
		rowDescVec[CouponRate] = couponRateDesc.str();
		rowTypeVec[CouponRate] = ARM_STRING;

		//FUNDING
		CC_Ostringstream fundDesc;
		fundDesc << "-FundRate[i]*FundLev[i]*Notional[i]";
		rowDescVec[Funding] = fundDesc.str();
		rowTypeVec[Funding] = ARM_STRING;

		//SUMCOUP
		CC_Ostringstream sumCoupDesc;
		ARM_Date fixingDate;
		double currFixing;
		double currCoupon;
		double coupons = 0.0;
		ARM_Date strikeDate;
		double currStrike;
		double currFixedRate;
		double currSpread;

		int sizePastFixings = 0;

		if (itsPastFixings)
		{
			sizePastFixings = itsPastFixings->size();
			for (int i = 0; i < sizePastFixings; i++)
			{
				fixingDate = itsPastFixings->GetDiscreteDates()->Elt(i);
				currFixing = itsPastFixings->GetDiscreteValues()->Elt(i);
				
				currStrike = const_cast< ARM_ReferenceValue& >(itsCapStrike).Interpolate(fixingDate);
				currFixedRate = const_cast< ARM_ReferenceValue& >(itsCapFixed).Interpolate(fixingDate);
				currSpread = const_cast< ARM_ReferenceValue& >(itsCapSpread).Interpolate(fixingDate);

				currCoupon = (currFixing < currStrike) ? currFixedRate : (currFixing+currSpread);
				
				coupons += currCoupon; 
			}
		}
		if (nbEvent > itsGlobalCapParams->Elt(0))
		{
			if (eventIdx < nbEvent-itsGlobalCapParams->Elt(0))
				sumCoupDesc << "0";
			else if (eventIdx == nbEvent-itsGlobalCapParams->Elt(0))
				sumCoupDesc << "CouponRate[i]+" << coupons;
			else
				sumCoupDesc << "SumCoup[i-1]+CouponRate[i]";
		}
		else
		{
			if(eventIdx == 0)
				sumCoupDesc << "CouponRate[i]+" << coupons;
			else
				sumCoupDesc << "SumCoup[i-1]+CouponRate[i]";
		}
		rowDescVec[SumCoup] = sumCoupDesc.str();
		rowTypeVec[SumCoup] = ARM_STRING;

		//COUPON
		CC_Ostringstream couponDesc;
		couponDesc << "CouponRate[i]*Basis[i]*DF[i]*Notional[i]";
		rowDescVec[Coupon] = couponDesc.str();
		rowTypeVec[Coupon] = ARM_STRING;

		//SWAP
		CC_Ostringstream swapDesc;
		swapDesc << "Funding[i]+Coupon[i]";
		rowDescVec[Swap] = swapDesc.str();
		rowTypeVec[Swap] = ARM_STRING;

		//GLOBALCAP
		CC_Ostringstream globalCapDesc;
		if(eventIdx == nbEvent-1)
		{
			globalCapDesc << "Max(SumCoup[i]/";
			globalCapDesc << itsGlobalCapParams->Elt(0) << "-";
			globalCapDesc << itsGlobalCapParams->Elt(2) << ",0)*";
			globalCapDesc << itsGlobalCapParams->Elt(1);
			globalCapDesc << "*Notional[i]*DF[i]";
		}
		else
			globalCapDesc << "0";
		rowDescVec[GlobalCap] = globalCapDesc.str();
		rowTypeVec[GlobalCap] = ARM_STRING;

		//PRODUCT
		CC_Ostringstream productDesc;
		productDesc << "Swap[i]-GlobalCap[i]";
		rowDescVec[Product] = productDesc.str();
		rowTypeVec[Product] = ARM_STRING;

		return ARM_RowInfo(rowDescVec,rowTypeVec);
	}
	catch (Exception& x)
	{
	    x.DebugPrint();

		throw x;
	}
	catch (...)
	{
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in GlobalCapCalculator::MiddleRows" );
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_GlobaCapCalculator
///	Routine: CreateAndSetModel
///	Returns: void
///	Action : create Markov Fonctional
/////////////////////////////////////////////////////////////////
void ARM_GlobalCapCalculator::CreateAndSetModel()
{
	/// Create MC Method

	// random generator
	ARM_RandomGenerator* randGen1=NULL;
	ARM_RandomGenerator* randGen2=NULL;
	ARM_RandomGenerator* randGen=NULL;

	ARM_RandomGeneratorPtrVector randGenVector;

	ARM_RandomGeneratorPtr pRandGen1, pRandGen2, pRandGen;
		
	randGen1 = ARM_RandGenFactory.Instance()->CreateRandGen(
			(ARM_RandGenFactoryImp::BaseGenType) ARM_ArgConv_BaseGenAlgoType.GetNumber(itsRandGenerator[0]));
	pRandGen1 = ARM_RandomGeneratorPtr( randGen1 );

	randGen2 = ARM_RandGenFactory.Instance()->CreateRandGen(
			ARM_RandGenFactoryImp::UnknownBaseGenAlgorithm,
			(ARM_RandGenFactoryImp::TransformAlgo) ARM_ArgConv_TransformAlgoType.GetNumber(itsRandGenerator[1]),
			pRandGen1);

	if (itsRandGenerator[2] != "NONE")
	{
		pRandGen2 = ARM_RandomGeneratorPtr( randGen2 );
		randGen = ARM_RandGenFactory.Instance()->CreateRandGen(
			ARM_RandGenFactoryImp::UnknownBaseGenAlgorithm,
			(ARM_RandGenFactoryImp::TransformAlgo)   ARM_ArgConv_TransformAlgoType.GetNumber(itsRandGenerator[2]),
			pRandGen2,
			 ARM_RandomGeneratorPtr(NULL),
			-1.0,
			1.0,
			1.0,
			ARM_GP_T_Vector<size_t>(1,10),
			4.0,
			0.0,
			(ARM_RandGenFactoryImp::RandGenOrder) ARM_ArgConv_RandGenOrder.GetNumber("PathOrder"));
	}
	else
		randGen = randGen2;

	pRandGen = ARM_RandomGeneratorPtr( randGen );

	if (itsRandGenerator.size() == 7)
	{
		ARM_RandomGenerator* randGen3=NULL;
		ARM_RandomGenerator* randGen4=NULL;
		ARM_RandomGenerator* randGen5=NULL;
		ARM_RandomGeneratorPtr pRandGen3, pRandGen4, pRandGen5;

		if (itsRandGenerator[5] != "NONE")
		{
			randGen3 = ARM_RandGenFactory.Instance()->CreateRandGen(
				ARM_RandGenFactoryImp::UnknownBaseGenAlgorithm,
				(ARM_RandGenFactoryImp::TransformAlgo)   ARM_ArgConv_TransformAlgoType.GetNumber(itsRandGenerator[5]),
				pRandGen);
			pRandGen3 = ARM_RandomGeneratorPtr( randGen3 );
		}
		else
			pRandGen3 = pRandGen;

		randGen4 = ARM_RandGenFactory.Instance()->CreateRandGen(
				(ARM_RandGenFactoryImp::BaseGenType) ARM_ArgConv_BaseGenAlgoType.GetNumber(itsRandGenerator[3]));
		pRandGen4 = ARM_RandomGeneratorPtr( randGen4 );

		randGen5 = ARM_RandGenFactory.Instance()->CreateRandGen(
				ARM_RandGenFactoryImp::UnknownBaseGenAlgorithm,
				(ARM_RandGenFactoryImp::TransformAlgo) ARM_ArgConv_TransformAlgoType.GetNumber(itsRandGenerator[4]),
				pRandGen4);
		pRandGen5 = ARM_RandomGeneratorPtr( randGen5 );

		randGen = ARM_RandGenFactory.Instance()->CreateRandGen(
				ARM_RandGenFactoryImp::UnknownBaseGenAlgorithm,
				(ARM_RandGenFactoryImp::TransformAlgo) ARM_ArgConv_TransformAlgoType.GetNumber("MixteGen"),
				pRandGen3,
				pRandGen5,
				-1,
				1,
				1,
				ARM_GP_T_Vector<size_t>(1,10),
				4.0,
				0,
				atoi(itsRandGenerator[6].c_str()));

		pRandGen = ARM_RandomGeneratorPtr( randGen );
	}

	randGenVector.resize(1);
	randGenVector[0] = pRandGen;

	// scheduler and sampler
	ARM_SchedulerBase* scheduler = NULL;
	ARM_SamplerBase* sampler = NULL;
	int multiDim = 2;
	std::vector<double> samplerDatas(0);

	int schedulerType = ARM_SchedulerBase::TimeStepPerYear;
	std::vector<double> schedulerDatas(1);
	schedulerDatas[0] = 1;

	int nbMCSteps;
	scheduler = ARM_SchedulerFactory.Instance()->CreateScheduler(
		schedulerType,
		schedulerDatas,
		nbMCSteps);

	sampler = ARM_SamplerFactory.Instance()->CreateSampler(
		multiDim,
		itsSamplerType,
		samplerDatas,
		scheduler);

	// mc method
	// will be set to model after calibration !
	itsMCMethod = ARM_NumMethodPtr(new ARM_MCMethod(itsNbSteps,
													randGenVector,
													sampler));

//	delete scheduler;
//	delete sampler;

	// model parameters
	ARM_CurveModelParam* meanReversion = dynamic_cast< ARM_CurveModelParam* >( GetMktDataManager()->GetData(GetKeys()[MrsKey]) );

	//int size = (*itsExerDateStrip).GetResetDates()->size();
	int size = (*itsExerDateStrip).GetResetDates()->size(); 
	int pastFixingsSize = 0;

	//if (itsPastFixings)
	//	pastFixingsSize = itsPastFixings->GetSize();
	
	//size -= pastFixingsSize;
	
	std::vector<double> breakPointTimes(size);
	std::vector<double> values(size);

	ARM_ZeroCurve* zc = dynamic_cast< ARM_BSSmiledModel* >( GetMktDataManager()->GetData(GetKeys()[CfModelKey]) )->GetZeroCurve();
	ARM_Date AsOf = zc->GetAsOfDate();

	for (int i=0; i<size; i++)
	{
		breakPointTimes[i] = (*itsExerDateStrip).GetResetDates()->Elt(i+pastFixingsSize) - AsOf.GetJulian();
		values[i] = itsCalibParams->Elt(3) + i*itsCalibParams->Elt(4);
	}

	ARM_CurveModelParam calibParamHKVol = ARM_CurveModelParam(ARM_ModelParamType::Volatility, 
																 &values, 
																 &breakPointTimes);

	vector<ARM_ModelParam*> paramsVector;
	paramsVector.push_back((ARM_ModelParam*)meanReversion);
	paramsVector.push_back((ARM_ModelParam*)&calibParamHKVol);

	// Create Model
	ARM_MarkovFunctional* mod= NULL;
	mod = new ARM_MarkovFunctional( CreateClonedPtr( zc ) );
	
	ARM_ModelParamVector modelParams(paramsVector);
    mod->SetModelParams(ARM_ModelParamsMF(modelParams));

	ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::TerminalZc ) );
	mod->SetNumeraire( numeraire );

	SetPricingModel(ARM_PricingModelPtr(static_cast< ARM_PricingModel* >(mod)));
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_GlobaCapCalculator
///	Routine: CreateAndSetCalibration
///	Returns: void
///	Action : create calibMethod
/////////////////////////////////////////////////////////////////
void ARM_GlobalCapCalculator::CreateAndSetCalibration()
{
	try
	{
		int size = (*itsExerDateStrip).GetResetDates()->size();

		ARM_StdPortfolio* CapletPf = new ARM_StdPortfolio();

		ARM_BSSmiledModel* capModel = dynamic_cast< ARM_BSSmiledModel* >( GetMktDataManager()->GetData(GetKeys()[CfModelKey]) );
		ARM_ZeroCurve* zc = capModel->GetZeroCurve();
		ARM_Date AsOf = zc->GetAsOfDate();
		int indexType = FromFrequencyToXiborType(itsPayFreq, zc->GetCurrencyUnit()->GetCcyName());

		string MethodTypeStr = "Numerical";
		ARM_MethodType methodType = (ARM_MethodType) ARM_ArgConv_CalibMethod.GetNumber(MethodTypeStr);

		string MktTargetStr = "UNKNOWN_TAR";
		ARM_MktTargetType  mktTargetType = (ARM_MktTargetType)ARM_ArgConv_TargetFuncMethod.GetNumber(MktTargetStr);

		ARM_VanillaSecDensityPtrVector vanillaSecDensities( size );

		ARM_VolCurve* alphaCurve = capModel->FromSigmaToAlpha();
		ARM_VolCurve* rhoCurve = capModel->GetRho();
		ARM_VolCurve* nuCurve = capModel->GetNu();
		ARM_VolCurve* betaCurve = capModel->GetBeta();

		// build another date strip for clibration schedule
		// when fixing is ARREARS, we should deal with FORWARD START and FORWARD END dates
		// instead of start and end dates !
		ARM_DateStrip* dateStrip = (ARM_DateStrip*)itsExerDateStrip->Clone();
		if (itsResetTiming == K_ARREARS)
		{
			dateStrip->SetFlowStartDates(dateStrip->GetFwdRateStartDates());
			dateStrip->SetFlowEndDates(dateStrip->GetFwdRateEndDates());
		}

		for (int i = 0; i < size; ++i)
		{
			ARM_VanillaSecurityDensity* vanillaSecDensity = NULL;
			// TMP : why don't we use DayCount ???
			double matu = ((*(dateStrip->GetResetDates()))[i] - AsOf.GetJulian()) / 365.0;
			if (matu < 0.0)
				matu = 0.0;
			double tenor = ((*(dateStrip->GetFlowEndDates()))[i] - (*(dateStrip->GetFlowStartDates()))[i]) / 365.0;
			
			double alpha = alphaCurve->ComputeVolatility(matu, tenor)/100.0;
			double rho = rhoCurve->ComputeVolatility(matu, tenor);
			double nu = nuCurve->ComputeVolatility(matu, tenor);
			double beta = betaCurve->ComputeVolatility(matu, tenor);

			ARM_SABRDensityFunctor* method = NULL;
			method = new ARM_SABRDensityFunctor(alpha, beta, rho, nu, capModel->GetSABRFlag());
			ARM_DensityFunctorPtr methodPtr(method);
			
			vanillaSecDensity = new ARM_VanillaSecurityDensity( (*(dateStrip->GetResetDates()))[i], 
																(*(dateStrip->GetFlowStartDates()))[i],
																(*(dateStrip->GetFlowEndDates()))[i],
																methodPtr ); 

			vanillaSecDensities[i] = ARM_VanillaSecDensityPtr( static_cast<ARM_VanillaSecurityDensity*> (vanillaSecDensity->Clone()) );
		}

		ARM_CalibMethod* calibMethod = new ARM_CalibMethod( ARM_StdPortfolioPtr(NULL),
															ARM_ModelParamVector(0),
															methodType, 
															100, 
															mktTargetType,
															NULL, 
															NULL, 
															false,
															0, /// C_factorNb
															1, /// NbIter
															true, 
															ARM_DateStripPtr(dateStrip),
															vanillaSecDensities);

		SetCalibMethod(ARM_CalibMethodPtr(calibMethod));
	}
	catch (Exception& x)
	{
	    x.DebugPrint();

		throw x;
	}
	catch (...)
	{
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in GlobalCapCalculator::CreateCalibMethods" );
	}
}

void ARM_GlobalCapCalculator::UpdateModel()
{
	CreateAndSetModel();
}

void ARM_GlobalCapCalculator::UpdateCalibration(bool isUpdateStrike)
{
	CreateAndSetCalibration();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_GlobaCapCalculator
///	Routine: Calibrate
///	Returns: void
///	Action : Calibrate the model using PDE Method
///          Set MC Method afterwards
/////////////////////////////////////////////////////////////////
void ARM_GlobalCapCalculator::Calibrate()
{
	try
	{
		// Create PDE Method
		ARM_PDEMethod* PDEMethod= NULL;

		// Gets the right NumScheme from its name
		int numSchemeType = ARM_ArgConv_PDENumSchemeType.GetNumber(string("CN1F"));
		ARM_PDENumericalScheme* numScheme = ARM_PDENumericalScheme::getNumericalSchemeInstanceById(numSchemeType);

		// Builds PDE NumMethod With this numscheme
		PDEMethod = new ARM_PDEMethod(numScheme,itsCalibParams->Elt(0), itsCalibParams->Elt(1), itsCalibParams->Elt(2));
		(&*GetPricingModel())->SetNumMethod(ARM_NumMethodPtr(PDEMethod));

		// Calibrate with PDE Method
		GetCalibMethod()->Calibrate(&*GetPricingModel());

		// Set MC Mtehod for pricing
		(&*GetPricingModel())->SetNumMethod(itsMCMethod);
	}
	catch (Exception& x)
	{
	    x.DebugPrint();

		throw x;
	}
	catch (...)
	{
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in GlobalCapCalculator::Calibrate" );
	}
}

ARM_RowInfo ARM_GlobalCapCalculator::ColumnNames() const
{
	// Number of Columns
	size_t colNamesSize;

	colNamesSize = sizeof(GlobalCapColNamesTable)/sizeof(GlobalCapColNamesTable[0]);

    vector< string > colNamesVec(colNamesSize);
    vector< ARM_GP_VALUE_TYPE > colTypeVec(colNamesSize, ARM_STRING); 

	for (size_t i = 0; i < colNamesSize; ++i)
	{
		colNamesVec[i] = GlobalCapColNamesTable[i];
	}

    ARM_RowInfo rowInfo(colNamesVec,colTypeVec);

    return rowInfo;
}

/*****************************************************************************
/ Class  : ARM_GlobalCapCalculator
/ Routine: PricedColumnNames
/ Returns: ARM_StringVector
/ Action : create the priced column names of the deal description
******************************************************************************/
ARM_StringVector ARM_GlobalCapCalculator::PricedColumnNames() const
{
	ARM_StringVector pricedColumns;

	// effective nb of columns to price :
	int colSize = MIN(NbProductsToPrice, itsProductsToPrice.size());

	for (int i=0; i<colSize; i++)
	{
		if (itsProductsToPrice[i])
			pricedColumns.push_back(GlobalCapColNamesTable[GlobalCapProductToPriceColumns[i]]);
	}

	return pricedColumns;
}

void ARM_GlobalCapCalculator::InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const
{
    string zeroValue("0");

	rowDescVec[Funding] = zeroValue;
	rowTypeVec[Funding] = ARM_DOUBLE;

	rowDescVec[GlobalCap] = zeroValue;
	rowTypeVec[GlobalCap] = ARM_DOUBLE;

	rowDescVec[Coupon] = zeroValue;
	rowTypeVec[Coupon] = ARM_DOUBLE;

	rowDescVec[Swap] = zeroValue;
	rowTypeVec[Swap] = ARM_DOUBLE;

	rowDescVec[Product] = zeroValue;
	rowTypeVec[Product] = ARM_DOUBLE;

	rowDescVec[SumCoup] = zeroValue;
	rowTypeVec[SumCoup] = ARM_DOUBLE;

	rowDescVec[DF] = zeroValue;
	rowTypeVec[DF] = ARM_DOUBLE;

	rowDescVec[FundRate] = zeroValue;
	rowTypeVec[FundRate] = ARM_DOUBLE;

	rowDescVec[Index] = zeroValue;
	rowTypeVec[Index] = ARM_DOUBLE;

}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_GlobalCapCalculator
///	Routine: Price
///	Returns: the global cap value
///	Action : price the global cap
/////////////////////////////////////////////////////////////////
double ARM_GlobalCapCalculator::Price()
{
	try
	{
		CalibrateAndTimeIt();

		/// Price the implicit product according to internal flag
		ARM_GenSecurityPtr genSec = GetGenSecurity();

		size_t nbPricedColumns = genSec->GetDealDescription().GetPricedColumnNames().size();

		double price = 0.0;
		ARM_GenPricer* genPricer = NULL;
		genPricer = new ARM_GenPricer( &*genSec, &*GetPricingModel() );
		ARM_AutoCleaner<ARM_GenPricer> HoldGP(genPricer );
		genPricer->Price();

		// effective nb of columns to price :
		int colSize = MIN(NbProductsToPrice, itsProductsToPrice.size());

		for (size_t i(0); i < colSize; ++i)
		{
			if (itsProductsToPrice[i])
			{
				price	= itsPayRec * genPricer->GetPricerInfo()->GetContents( GlobalCapColNamesTable[ GlobalCapProductToPriceColumns[i] ] ).GetData("Price").GetDouble();
						
				if (i == FundingPrice)
				{
					itsFundingPrice = price;
				}
				else if (i == GlobalCapPrice)
				{
					itsGlobalCapPrice = price;
				}
				else if (i == CouponPrice)
				{
					itsCouponPrice = price;
				}
				else if (i == SwapPrice)
				{
					itsSwapPrice = price;
				}
				else if (i == ProductPrice)
				{
					itsProductPrice = price;
				}
				else if (i == SumCoupPrice)
				{
					itsSumCoupPrice = price;
				}
				else if (i == DFPrice)
				{
					itsDFPrice = price;
				}
				else if (i == FundRatePrice)
				{
					itsFundRatePrice = price;
				}
				else if (i == IndexPrice)
				{
					itsIndexPrice = price;
				}
			}
			itsHasBeenPriced = true;
		}

		itsProductsToPrice[0] ? price = itsProductPrice : itsGlobalCapPrice ;
		return price;
	}
	catch (Exception& x)
	{
	    x.DebugPrint();

		throw x;
	}
	catch (...)
	{
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in GlobalCapCalculator::Price" );
	}
}

/*****************************************************************************
/ Class   : ARM_GlobalCapCalculator
/ Routines: ComputePricingData
/ Returns :
/ Action  : set pricing results into itsPricingData
******************************************************************************/
void ARM_GlobalCapCalculator::ComputePricingData() const
{
	if (!itsHasBeenPriced)
		const_cast<ARM_GlobalCapCalculator*>(this)->PriceAndTimeIt();

	GetPricingData()[ "Funding" ]			= itsFundingPrice;
	GetPricingData()[ "Coupon" ]			= itsCouponPrice;
	GetPricingData()[ "Swap" ]				= itsSwapPrice;		
	GetPricingData()[ "Product" ]			= itsProductPrice;		
	GetPricingData()[ "GlobalCap" ]			= itsGlobalCapPrice;	
	GetPricingData()[ "SumCoup" ]			= itsSumCoupPrice;	
	GetPricingData()[ "DF" ]				= itsDFPrice;	
	GetPricingData()[ "FundRate" ]			= itsFundRatePrice;	
	GetPricingData()[ "Index" ]				= itsIndexPrice;	
}

/*****************************************************************************
/ Class   : ARM_GlobalCapCalculator
/ Routines: GeneralDataToString
/ Returns :
/ Action  : Construct information string for
/           general product data
******************************************************************************/
string ARM_GlobalCapCalculator::GeneralDataToString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream generalGCData;

	//General Data
	generalGCData << "\n\nGENERAL DATA"	<< endl;
	generalGCData << "Currency\t\t: "	<<  itsCcy.GetCcyName() << endl;
	generalGCData << "StartDate\t\t: "  <<  itsStartDate.toString() << endl;
    generalGCData << "EndDate\t\t\t: "	<<  itsEndDate.toString() << endl;
    generalGCData << "Pay/Rec\t\t\t: "	<<  ARM_ParamView::GetMappingName(S_RECEIVE_PAY, itsPayRec) << endl;
	generalGCData << "\n";
	
	//Fund Data
	generalGCData << "FUND DATA" << endl;
	generalGCData << "Day Count\t\t: "	<< ARM_ArgConvReverse_DayCount.GetString(itsDayCount) << "\n";
	generalGCData << "Frequency\t\t: "	<< ARM_ArgConvReverse_StdFrequency.GetString(itsPayFreq) << endl;
   	generalGCData << "Reset Gap\t\t: "	<< itsResetGap << "\n";
	generalGCData << "Pay Gap\t\t: "	<< itsPayGap << "\n";
	generalGCData << "Pay Timing\t\t: "	<< itsPayTiming << "\n";
	generalGCData << "\n";

	generalGCData << endl;

    return generalGCData.str();
}

/*****************************************************************************
/ Class   : ARM_GlobalCapCalculator
/ Routines: DealDesDataToString
/ Returns :
/ Action  : Construct information string for deal description data
******************************************************************************/
string ARM_GlobalCapCalculator::DealDesDataToString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream dealDesGCData;

	double asOf = GetMktDataManager()->GetAsOfDate().GetJulian();
	//Notional Profile
	dealDesGCData << "=====================Notional Profile===================" << endl;
	ARM_GP_CurvePtr notional(RefValueToCurve(itsNotional, asOf));
	dealDesGCData << notional->toString();

	//Fund Lev Profile
	dealDesGCData << "\n=====================Fund Lev Profile===================" << endl;
	ARM_GP_CurvePtr fundLev(RefValueToCurve(itsFundLev, asOf));
	dealDesGCData << fundLev->toString();

	//Cap Lev Profile
	dealDesGCData << "\n======================Cap Lev Profile===================" << endl;
	ARM_GP_CurvePtr capLev(RefValueToCurve(itsCapLev, asOf));
	dealDesGCData << capLev->toString();

	//Cap Fixed Profile
	dealDesGCData << "\n=====================Cap Fixed Profile==================" << endl;
	ARM_GP_CurvePtr capFixed(RefValueToCurve(itsCapSpread, asOf));
	dealDesGCData << capFixed->toString();

	//Cap Strike Profile
	dealDesGCData << "\n=====================Cap Strike Profile=================" << endl;
	ARM_GP_CurvePtr capStrike(RefValueToCurve(itsCapStrike, asOf));
	dealDesGCData << capStrike->toString();

	//Cap Spread Profile
	dealDesGCData << "\n=====================Cap Spread Profile=================" << endl;
	ARM_GP_CurvePtr capSpread(RefValueToCurve(itsCapSpread, asOf));
	dealDesGCData << capSpread->toString();

    return dealDesGCData.str();
}

/*****************************************************************************
/ Class  : ARM_GlobalCapCalculator
/ Routine: View
/ Returns: 
/ Action : .
******************************************************************************/
void ARM_GlobalCapCalculator::View(char* id, FILE* ficOut) const
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

    fprintf(fOut,"\n========================================================");
	fprintf(fOut,"\n===============> GLOBAL CAP CALCULATOR <================");
    fprintf(fOut,"\n========================================================");

	fprintf(fOut,"%s",GeneralDataToString().c_str());
	fprintf(fOut,"%s",DealDesDataToString().c_str());

	CC_Ostringstream genData;
	fprintf(fOut,"\n\n ============>    GEN CALCULATOR    <=========== \n");
	genData << "\n\nCommon Calculator Part\n" << ARM_GenCalculator::toString() << "\n\n";
	fprintf(fOut,"%s",genData.str().c_str());
	fprintf(fOut,"\n\n ============>  EOF GEN CALCULATOR  <=========== \n");

	if ( ficOut == NULL )
       fclose(fOut);
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
