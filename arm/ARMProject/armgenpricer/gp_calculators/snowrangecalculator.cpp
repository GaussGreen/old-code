/*!
 /** Copyright (c) CDC IXIS CM July 2003 Paris
 *	\file SnowRangeCalculator.cpp
 *  \brief file for the Snow Range Calculator
 *	\author  P. LAM
 *	\version 1.0
 *	\date Mar 2006
 */

/// this header comes first as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpcalculators/snowrangecalculator.h"
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
#include "gpmodels/SmiledFRM.h"

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
#include "gpnummethods/amc_ls.h"

/// kernel
#include "inst/swapleg.h"
#include "inst/fixleg.h"
#include "inst/portfolio.h"
#include "inst/barrier.h"

CC_BEGIN_NAMESPACE( ARM )

const int NB_FACTORS				= 1;

const unsigned int SNOWRANGE_SCHED	= 0;

const int DEAL_NBCOLUMNS			= 28;

const string ARM_SnowRangeCalculator::SnowRangeColNamesTable [] =
{  
	"ResetDate",
	"StartDate",
	"EndDate",
	"DealEndDate",
	"PayDate",
	"Notional",
	"Spread",
	"Strike",
	"Ratchet",
	"CF", 
	"FundingBasis",
	"CouponBasis",
	"Num",
	"Index",
	"Libor",
	"DF",
	"Funding",
	"RatioADV",
	"RatioARR",
	"NbPeriod",
	"FCorr",
	"CouponLev",
	"Coupon",
	"StdCoupon",
	"PlainCoupon",
	"Swap",
	"CallOption",
	"CallSwap",
};


const int ARM_SnowRangeCalculator::SnowRangeProductToPriceColumns [] =
{
	CallSwap,
	CallOption,
	Swap,
	PlainCoupon,
	StdCoupon,
	Coupon,
	Funding,
	Libor,
	Index,
	DF,
	RatioADV,
	RatioARR,
	FCorr,
	CouponLev,
};


/*****************************************************************************
/ Class  : ARM_SnowRangeCalculator
/ Routines: ARM_SnowRangeCalculator
/ Returns :
/ Action  : basic constructor
******************************************************************************/
ARM_SnowRangeCalculator::ARM_SnowRangeCalculator( ARM_Currency& ccy,
												  ARM_Date&	startDate,
												  ARM_Date&	endDate,
												  int payRec,
												  ARM_ReferenceValue& notional,
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
												  ARM_ReferenceValue& spread,
												  ARM_ReferenceValue& strike,
												  ARM_ReferenceValue& ratchet,
												  ARM_ReferenceValue& cashFlow,
												  ARM_ReferenceValue& fixedRate,
												  ARM_ReferenceValue& leverage,
												  ARM_Vector* snowRangeParams,
												  ARM_Vector* calibParams,
												  string modelName,
												  ARM_Vector* modelParams,
												  int nbSteps,
												  string generatorType,
												  string inversionMethod,
												  bool antithetic,
												  int samplerType,
												  ARM_StringVector& mdmKeys,
												  const ARM_MarketData_ManagerRep& mktDataManager,
												  const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice)
: ARM_GenCalculator(mktDataManager),
itsCcy(ccy),
itsStartDate(startDate),
itsEndDate(endDate),
itsSchedEndDate(endDate),
itsPayRec(payRec),
itsNotional(notional),
itsFundingIndexTerm(fundingIndexTerm),
itsFundingDayCount(fundingDayCount),
itsCouponIndexTerm(couponIndexTerm),
itsCouponDayCount(couponDayCount),
itsResetFreq(resetFreq),
itsPayFreq(payFreq),
itsResetTiming(resetTiming),
itsPayTiming(payTiming),
itsResetGap(resetGap),
itsResetCal(resetCal),
itsPayCal(payCal),
itsAdjRule(adjRule),
itsIntRule(intRule),
itsSpread(spread),
itsStrike(strike),
itsRatchet(ratchet),
itsCashFlow(cashFlow),
itsFixedRate(fixedRate),
itsLeverage(leverage),
itsSnowRangeParams(snowRangeParams),
itsCalibParams(calibParams),
itsModelName(modelName),
itsModelParams(modelParams),
itsNbSteps(nbSteps),
itsGeneratorType(generatorType),
itsInversionMethod(inversionMethod),
itsAntithetic(antithetic),
itsSamplerType(samplerType),
itsProductsToPrice(productsToPrice),
itsHasBeenPriced(false),
itsProductDateStrip(ARM_DateStripPtr(NULL)),
itsCalibDateStrip(ARM_DateStripPtr(NULL)),
itsCallSwapPrice(0.0),
itsCallOptionPrice(0.0),
itsSwapPrice(0.0),
itsPlainCouponPrice(0.0),
itsStdCouponPrice(0.0),
itsCouponPrice(0.0),
itsFundingPrice(0.0),
itsLiborPrice(0.0),
itsIndexPrice(0.0),
itsDFPrice(0.0),
itsRatioADVPrice(0.0),
itsRatioARRPrice(0.0),
itsFCorrPrice(0.0),	
itsCouponLevPrice(0.0)
{
	// WARNING :
	// Index rate is fixed for the next period.
	// 1st period : fixing but no payment is made for this first period.
	if (itsResetTiming == K_ARREARS)
		itsSchedEndDate.AddPeriod(itsCouponIndexTerm, const_cast<char*>(itsResetCal.c_str())).GetJulian();

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
/ Class  : ARM_SnowRangeCalculator
/ Routines: 
/ Returns :
/ Action  : copy constructor
******************************************************************************/
ARM_SnowRangeCalculator::ARM_SnowRangeCalculator(const ARM_SnowRangeCalculator& rhs)
: ARM_GenCalculator(rhs)
{
	CopyNoCleanUp(rhs);
}
	
/*****************************************************************************
/ Class  : ARM_SnowRangeCalculator
/ Routine: CopyNoCleanUp
/ Returns: void
/ Action : Arguments copy
*****************************************************************************/
void ARM_SnowRangeCalculator::CopyNoCleanUp(const ARM_SnowRangeCalculator& rhs)
{
	itsCcy							= rhs.itsCcy;
	itsStartDate					= rhs.itsStartDate;        
	itsEndDate						= rhs.itsEndDate;    
	itsSchedEndDate					= rhs.itsSchedEndDate;    
	itsPayRec						= rhs.itsPayRec;
	itsNotional						= rhs.itsNotional;
	itsFundingDayCount				= rhs.itsFundingDayCount;
	itsCouponDayCount				= rhs.itsCouponDayCount;
	itsResetFreq					= rhs.itsResetFreq;
	itsPayFreq						= rhs.itsPayFreq;
	itsResetGap						= rhs.itsResetGap;
	itsPayGap						= rhs.itsPayGap;
	itsResetTiming					= rhs.itsResetTiming;
	itsPayTiming					= rhs.itsPayTiming;
	itsResetCal						= rhs.itsResetCal;
	itsPayCal						= rhs.itsPayCal;
	itsAdjRule						= rhs.itsAdjRule;
	itsIntRule						= rhs.itsIntRule;
	itsSpread						= rhs.itsSpread;
	itsFixedRate					= rhs.itsFixedRate;
	itsLeverage						= rhs.itsLeverage;
	itsStrike						= rhs.itsStrike;
	itsRatchet						= rhs.itsRatchet;
	itsCashFlow						= rhs.itsCashFlow;
	itsSnowRangeParams				= rhs.itsSnowRangeParams ? (ARM_Vector*)rhs.itsSnowRangeParams->Clone() : NULL;
	itsCalibParams					= rhs.itsCalibParams ? (ARM_Vector*)rhs.itsCalibParams->Clone() : NULL;
	itsModelName					= rhs.itsModelName;
	itsModelParams					= rhs.itsModelParams ? (ARM_Vector*)rhs.itsModelParams->Clone() : NULL;
	itsNbSteps						= rhs.itsNbSteps;
	itsGeneratorType				= rhs.itsGeneratorType;
	itsInversionMethod				= rhs.itsInversionMethod;
	itsAntithetic					= rhs.itsAntithetic;
	itsSamplerType					= rhs.itsSamplerType;
	itsProductsToPrice				= rhs.itsProductsToPrice;
	itsHasBeenPriced				= rhs.itsHasBeenPriced;
	itsCallSwapPrice				= rhs.itsCallSwapPrice;
	itsCallOptionPrice				= rhs.itsCallOptionPrice;
	itsSwapPrice					= rhs.itsSwapPrice;
	itsPlainCouponPrice				= rhs.itsPlainCouponPrice;
	itsStdCouponPrice				= rhs.itsStdCouponPrice;
	itsCouponPrice					= rhs.itsCouponPrice;
	itsFundingPrice					= rhs.itsFundingPrice;
	itsLiborPrice					= rhs.itsLiborPrice;
	itsIndexPrice					= rhs.itsIndexPrice;
	itsDFPrice						= rhs.itsDFPrice;
	itsProductDateStrip				= (ARM_DateStripPtr((!rhs.itsProductDateStrip.IsNull()?static_cast<ARM_DateStrip*>(rhs.itsProductDateStrip->Clone()):NULL)));
	itsCalibDateStrip				= (ARM_DateStripPtr((!rhs.itsCalibDateStrip.IsNull()?static_cast<ARM_DateStrip*>(rhs.itsCalibDateStrip->Clone()):NULL)));
}

/*****************************************************************************
/ Class  : ARM_SnowRangeCalculator
/ Routines: Clone
/ Returns :
/ Action  : Standard ARM object support
/           Call copy constructor
******************************************************************************/
ARM_Object*	ARM_SnowRangeCalculator::Clone() const
{
	return(new ARM_SnowRangeCalculator(*this));
}
	
/*****************************************************************************
/ Class  : ARM_SnowRangeCalculator
/ Routines: operator=
/ Returns :
/ Action  : Call copy constructor
******************************************************************************/
ARM_SnowRangeCalculator& ARM_SnowRangeCalculator::operator=(const ARM_SnowRangeCalculator& rhs)
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
/ Class  : ARM_SnowRangeCalculator
/ Routine: CleanUp
/ Returns: void
/ Action : delete pointors
******************************************************************************/
void ARM_SnowRangeCalculator::CleanUp()
{
	if (itsSnowRangeParams)
		delete itsSnowRangeParams;
	itsSnowRangeParams = NULL;

	if (itsCalibParams)
		delete itsCalibParams;
	itsCalibParams = NULL;

	if (itsModelParams)
		delete itsModelParams;
	itsModelParams = NULL;
}

/*****************************************************************************
/ Class  : ARM_SnowRangeCalculator
/ Routine: Destructor
/ Returns: void
/ Action : destroys the object
******************************************************************************/
ARM_SnowRangeCalculator::~ARM_SnowRangeCalculator()
{
	CleanUp();
}

/*****************************************************************************
/ Class  : ARM_SnowRangeCalculator
/ Routine: CheckMktData
/ Returns: 
/ Action : Checks if market data are of good type
******************************************************************************/
void ARM_SnowRangeCalculator::CheckMktData()
{
    /// MdM datas checking
	ARM_ZeroCurve* ycCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));
    if(!ycCurve)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Yield Curve for key=" + GetKeys()[YcKey] + " is expected in the Market Data Manager");

	ARM_BSModel* cfBSModel = dynamic_cast< ARM_BSModel* >(GetMktDataManager()->GetData(GetKeys()[CfModelKey]) );
    if(!cfBSModel)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : cap/floor B&S smiled model for key=" + GetKeys()[CfModelKey] + " is expected in the Market Data Manager");

	ARM_BSModel* swoptBSModel = dynamic_cast< ARM_BSModel* >(GetMktDataManager()->GetData(GetKeys()[OswModelKey]) );
    if(!swoptBSModel)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : swaption B&S smiled model for key=" + GetKeys()[OswModelKey] + " is expected in the Market Data Manager");

	ARM_ModelParam* mrsParam = dynamic_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[MrsKey]));
    if(!mrsParam || mrsParam->GetType() != ARM_ModelParamType::MeanReversion)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : MRS Param for key=" + GetKeys()[MrsKey] + " is expected in the Market Data Manager");

	ARM_ModelParam* humpParam = dynamic_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[HumpKey]));
    if(!humpParam || humpParam->GetType() != ARM_ModelParamType::Hump)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Hump Param for key=" + GetKeys()[HumpKey] + " is expected in the Market Data Manager");
}

/*****************************************************************************
/ Class  : ARM_SnowRangeCalculator
/ Routine: CheckData
/ Returns: 
/ Action : For the moment, only checks market data
******************************************************************************/
void ARM_SnowRangeCalculator::CheckData()
{
	CheckMktData();
}

/*****************************************************************************
/ Class  : ARM_SnowRangeCalculator
/ Routine: CreateCstManager
/ Returns: void
/ Action : create the const manager (static data).
******************************************************************************/
ARM_CstManagerPtr ARM_SnowRangeCalculator::CreateCstManager()
{
	try
	{
		double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

		vector<string>	cstNames; 

		cstNames.push_back("NotionalCst");
		cstNames.push_back("FixedRateCst");
		cstNames.push_back("LeverageCst");
		cstNames.push_back("StrikeCst");
		cstNames.push_back("SpreadCst");
		cstNames.push_back("RatchetCst");
		cstNames.push_back("CashFlowCst");

		vector<ARM_GramFctorArg> cstVector;
		
		//Notional
		ARM_Curve* tmpNotional = RefValueToCurve(itsNotional, asOfDate);
		cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpNotional))));

		//Fixed Rate
		ARM_Curve* tmpFixedRate = RefValueToCurve(itsFixedRate, asOfDate);
		cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpFixedRate))));

		//Leverage
		ARM_Curve* tmpLeverage = RefValueToCurve(itsLeverage, asOfDate);
		cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpLeverage))));

		//Strike
		ARM_Curve* tmpStrike = RefValueToCurve(itsStrike, asOfDate);
		cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpStrike))));

		//Spread
		ARM_Curve* tmpSpread = RefValueToCurve(itsSpread, asOfDate);
		cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpSpread))));

		//Ratchet
		ARM_Curve* tmpRatchet = RefValueToCurve(itsRatchet, asOfDate);
		cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpRatchet))));

		//CF
		ARM_Curve* tmpCashFlow = RefValueToCurve(itsCashFlow, asOfDate);
		cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpCashFlow))));

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
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in SnowRangeCalculator::CreateCstManager" );
	}
}

	//Outputs	
ARM_RowInfo	ARM_SnowRangeCalculator::ColumnNames() const
{
	// Number of Columns
	size_t colNamesSize;

	colNamesSize = sizeof(SnowRangeColNamesTable)/sizeof(SnowRangeColNamesTable[0]);

    vector< string > colNamesVec(colNamesSize);
    vector< ARM_GP_VALUE_TYPE > colTypeVec(colNamesSize, ARM_STRING); 

	for (size_t i = 0; i < colNamesSize; ++i)
	{
		colNamesVec[i] = SnowRangeColNamesTable[i];
	}

    ARM_RowInfo rowInfo(colNamesVec,colTypeVec);

    return rowInfo;
}

	// Some utilities
//void ARM_SnowRangeCalculator::GenerateProductDescription(ARM_Portfolio* portfolio)
//{}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_SnowRangeCalculator
///	Routine: CustomDatesStructure
///	Returns: ARM_DateStripCombiner
///	Action : create the list of all event dates of the Snow Range.
///		customized Schedule Not implemented yet.
/////////////////////////////////////////////////////////////////
ARM_DateStripCombiner ARM_SnowRangeCalculator::CustomDatesStructure(const ARM_DateStripVector& dateStrips) const
{
	return ARM_DateStripCombiner();
}


/*****************************************************************************
/ Class  : ARM_SnowRangeCalculator
/ Routine: DatesStructure
/ Returns: ARM_DateStripVector
/ Action : create the list of all event dates of the Snow Range. 
******************************************************************************/
ARM_DateStripCombiner ARM_SnowRangeCalculator::DatesStructure() const
{
	try
	{
		ARM_DateStripVector SchedVect(1,NULL);
	
		double startDate = itsStartDate.GetJulian();
		double endDate = itsSchedEndDate.GetJulian();

		// calib date Strip
		ARM_DateStrip CalibSched(startDate, endDate, itsResetFreq, itsCouponDayCount,
								 itsResetCal.c_str(), itsAdjRule, itsIntRule, 
								 K_SHORTSTART, itsResetGap, itsResetFreq, GETDEFAULTVALUE,
								 itsPayCal.c_str(), K_ADVANCE, itsPayTiming, true);

		itsCalibDateStrip = ARM_DateStripPtr(static_cast<ARM_DateStrip*>(CalibSched.Clone()));
	
		// product date strip
		// last reset date can not be beyond product's end date
		ARM_Date tmpEndDate(itsEndDate);
		if (itsSchedEndDate > tmpEndDate.AddPeriod(itsResetFreq, const_cast<char*>(itsResetCal.c_str())))
		{
			endDate = tmpEndDate.GetJulian();
		}
		// if (IndexTerm > Freq), we have to adjust the fwd end dates
		ARM_DateStrip ProdSched(startDate, endDate, itsResetFreq, itsCouponDayCount,
								itsResetCal.c_str(), itsAdjRule, itsIntRule, 
								K_SHORTSTART, itsResetGap, itsResetFreq, GETDEFAULTVALUE,
								itsPayCal.c_str(), K_ADVANCE, itsPayTiming, true);

		if (StringMatuToYearTerm(const_cast<char*>(itsCouponIndexTerm.c_str())) > 1.0/itsResetFreq)
		{
			for (int i=0; i<ProdSched.GetResetDates()->size(); i++)
			{
				ARM_Date tmpStartDate(ProdSched.GetFwdRateStartDates()->Elt(i));			
				tmpStartDate.AddPeriod(itsCouponIndexTerm, const_cast<char*>(itsResetCal.c_str()));
				ProdSched.GetFwdRateEndDates()->Elt(i) = tmpStartDate.AdjustToBusDate(const_cast<char*>(itsResetCal.c_str()), itsIntRule).GetJulian();
			}
		}

		itsProductDateStrip = ARM_DateStripPtr(static_cast<ARM_DateStrip*>(ProdSched.Clone()));

		SchedVect[SNOWRANGE_SCHED] = &ProdSched;

		return SchedVect;
	}
	catch (Exception& x)
	{
	    x.DebugPrint();

		throw x;
	}
	catch (...)
	{
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in SnowRangeCalculator::DatesStructure" );
	}
}


/*****************************************************************************
/ Class  : ARM_SnowRangeCalculator
/ Routine: MiddleRows
/ Returns: ARM_RowInfo
/ Action : create a row of a deal description
/          warning : no payoff for the first period !!
*****************************************************************************/
ARM_RowInfo ARM_SnowRangeCalculator::MiddleRows(size_t eventIdx, 
                                          const ARM_DateStripCombiner& datesStructure) const
{
	try
	{
		//Useful data
		char* ccy	= itsCcy.GetCcyName();

		//Number of columns
		size_t descSize = DEAL_NBCOLUMNS;
	
		size_t nbEvent = datesStructure.GetDateStrip(SNOWRANGE_SCHED)->GetResetDates()->size();
		vector< string > rowDescVec(descSize);
		vector< ARM_GP_VALUE_TYPE > rowTypeVec(descSize, ARM_MISSING_TYPE); 
		
		InitPriceableColumns(rowDescVec,rowTypeVec);

		//RESET DATE
		double resetDate = (*(datesStructure.GetDateStrip(SNOWRANGE_SCHED)->GetResetDates()))[eventIdx];
		CC_Ostringstream resetDateDesc;
		resetDateDesc << CC_NS(std,fixed) << resetDate;
		rowDescVec[ResetDate] = resetDateDesc.str();
		rowTypeVec[ResetDate] = ARM_DATE_TYPE;

		//START DATE (in fact : FWD START ! so that we take Arrears fixings into account)
		double startDate = (*(datesStructure.GetDateStrip(SNOWRANGE_SCHED)->GetFwdRateStartDates()))[eventIdx] ;
		CC_Ostringstream startDateDesc;
		startDateDesc << CC_NS(std,fixed) << startDate;
		rowDescVec[StartDate] = startDateDesc.str();
		rowTypeVec[StartDate] = ARM_DATE_TYPE;

		//END DATE (in fact : FWD END ! so that we take Arrears fixings into account)
		// if the calculated end date is beyong product's real end date, take this real end date !
		double endDate = (*(datesStructure.GetDateStrip(SNOWRANGE_SCHED)->GetFwdRateEndDates()))[eventIdx];
		if (ARM_Date(endDate) > itsSchedEndDate)
			endDate = itsSchedEndDate.GetJulian();

		CC_Ostringstream endDateDesc;
		endDateDesc << CC_NS(std,fixed) << endDate;
		rowDescVec[EndDate] = endDateDesc.str();
		rowTypeVec[EndDate] = ARM_DATE_TYPE;

		//DEAL END DATE (for numeraire in AMC regression)
		CC_Ostringstream dealEndDateDesc;
		dealEndDateDesc << CC_NS(std,fixed) << itsSchedEndDate.GetJulian();
		rowDescVec[DealEndDate] = dealEndDateDesc.str();
		rowTypeVec[DealEndDate] = ARM_DATE_TYPE;

		//PAY DATE
		CC_Ostringstream payDateDesc;
		if (eventIdx == 0)
		{
			payDateDesc << CC_NS(std,fixed) << "na";
			rowDescVec[PayDate] = payDateDesc.str();
			rowTypeVec[PayDate] = ARM_STRING;
		}
		else
		{
			double payDate = (*(datesStructure.GetDateStrip(SNOWRANGE_SCHED)->GetPaymentDates()))[eventIdx-1] ;
			payDateDesc << CC_NS(std,fixed) << payDate;
			rowDescVec[PayDate] = payDateDesc.str();
			rowTypeVec[PayDate] = ARM_DATE_TYPE;
		}
		
		//NOTIONAL
		CC_Ostringstream notionalDesc;
		if (eventIdx == 0)
		{
			notionalDesc << CC_NS(std,fixed) << "na";
			rowDescVec[Notional] = notionalDesc.str();
			rowTypeVec[Notional] = ARM_STRING;
		}
		else
		{
			double notional = const_cast< ARM_ReferenceValue& >(itsNotional).Interpolate(startDate);
			notionalDesc << CC_NS(std,fixed) << notional;
			rowDescVec[Notional] = notionalDesc.str();
			rowTypeVec[Notional] = ARM_DOUBLE;
		}

		//SPREAD
		CC_Ostringstream spreadDesc;
		if (eventIdx == 0)
		{
			spreadDesc << CC_NS(std,fixed) << "na";
			rowDescVec[Spread] = spreadDesc.str();
			rowTypeVec[Spread] = ARM_STRING;
		}
		else
		{
			double spread = const_cast< ARM_ReferenceValue& >(itsSpread).Interpolate(startDate);
			spreadDesc << CC_NS(std,fixed) << spread;
			rowDescVec[Spread] = spreadDesc.str();
			rowTypeVec[Spread] = ARM_DOUBLE;
		}

		//STRIKE
		CC_Ostringstream strikeDesc;
		if (eventIdx == 0)
		{
			strikeDesc << CC_NS(std,fixed) << "na";
			rowDescVec[Strike] = strikeDesc.str();
			rowTypeVec[Strike] = ARM_STRING;
		}
		else
		{
			double strike = const_cast< ARM_ReferenceValue& >(itsStrike).Interpolate(startDate);
			strikeDesc << CC_NS(std,fixed) << strike;
			rowDescVec[Strike] = strikeDesc.str();
			rowTypeVec[Strike] = ARM_DOUBLE;
		}

		//RATCHET
		CC_Ostringstream ratchetDesc;
		if (eventIdx == 0)
		{
			ratchetDesc << CC_NS(std,fixed) << "na";
			rowDescVec[Ratchet] = ratchetDesc.str();
			rowTypeVec[Ratchet] = ARM_STRING;
		}
		else
		{
			double ratchet = const_cast< ARM_ReferenceValue& >(itsRatchet).Interpolate(startDate);
			ratchetDesc << CC_NS(std,fixed) << ratchet;
			rowDescVec[Ratchet] = ratchetDesc.str();
			rowTypeVec[Ratchet] = ARM_DOUBLE;
		}

		//CASHFLOW
		CC_Ostringstream cashFlowDesc;
		if (eventIdx == 0)
		{
			cashFlowDesc << CC_NS(std,fixed) << "na";
			rowDescVec[CashFlow] = cashFlowDesc.str();
			rowTypeVec[CashFlow] = ARM_STRING;
		}
		else
		{
			double cashFlow = const_cast< ARM_ReferenceValue& >(itsCashFlow).Interpolate(startDate);
			cashFlowDesc << CC_NS(std,fixed) << cashFlow;
			rowDescVec[CashFlow] = cashFlowDesc.str();
			rowTypeVec[CashFlow] = ARM_DOUBLE;
		}

		//FUNDING BASIS
		CC_Ostringstream fundingBasisDesc;
		if (eventIdx == 0)
		{
			fundingBasisDesc << CC_NS(std,fixed) << "na";
			rowDescVec[FundingBasis] = fundingBasisDesc.str();
			rowTypeVec[FundingBasis] = ARM_STRING;
		}
		else
		{
			// PREVIOUS START and PAY DATES !!!!!! (for basis computation)
			double prevStart = (*(datesStructure.GetDateStrip(SNOWRANGE_SCHED)->GetFwdRateStartDates()))[eventIdx-1] ;
			double prevPay = (*(datesStructure.GetDateStrip(SNOWRANGE_SCHED)->GetPaymentDates()))[eventIdx-1] ;

			double fundingBasis = CountYears(itsFundingDayCount, prevStart, prevPay);
			fundingBasisDesc << CC_NS(std,fixed) << fundingBasis;
			rowDescVec[FundingBasis] = fundingBasisDesc.str();
			rowTypeVec[FundingBasis] = ARM_DOUBLE;
		}

		//COUPON BASIS
		CC_Ostringstream couponBasisDesc;
		if (eventIdx == 0)
		{
			couponBasisDesc << CC_NS(std,fixed) << "na";
			rowDescVec[CouponBasis] = couponBasisDesc.str();
			rowTypeVec[CouponBasis] = ARM_STRING;
		}
		else
		{
			// PREVIOUS START and PAY DATES !!!!!! (for basis computation)
			double prevStart = (*(datesStructure.GetDateStrip(SNOWRANGE_SCHED)->GetFwdRateStartDates()))[eventIdx-1] ;
			double prevPay = (*(datesStructure.GetDateStrip(SNOWRANGE_SCHED)->GetPaymentDates()))[eventIdx-1] ;

			double couponBasis = CountYears(itsCouponDayCount, prevStart, prevPay);
			couponBasisDesc << CC_NS(std,fixed) << couponBasis;
			rowDescVec[CouponBasis] = couponBasisDesc.str();
			rowTypeVec[CouponBasis] = ARM_DOUBLE;
		}

		//NUMERAIRE
		CC_Ostringstream numDesc;
		if (eventIdx == 0)
		{
			numDesc << "0.0";
		}
		else
		{
			numDesc << "DF(" << ccy << ",";
			numDesc << "DealEndDate[i])";
		}
		rowDescVec[Numeraire] = numDesc.str();
		rowTypeVec[Numeraire] = ARM_STRING;
		
		//INDEX
		CC_Ostringstream indexDesc;
		indexDesc << "LIBOR(" << ccy << ",";
		indexDesc << "StartDate[i],";
		indexDesc << itsCouponIndexTerm << ")";
		rowDescVec[Index] = indexDesc.str();
		rowTypeVec[Index] = ARM_STRING;

		//LIBOR
		CC_Ostringstream liborDesc;
		liborDesc << "LIBOR(" << ccy << ",";
		liborDesc << "StartDate[i],";
		liborDesc << itsFundingIndexTerm << ")";
		rowDescVec[Libor] = liborDesc.str();
		rowTypeVec[Libor] = ARM_STRING;

		//DF
		CC_Ostringstream dfDesc;
		if (eventIdx == 0)
		{
			dfDesc << "0.0";
		}
		else
		{
			dfDesc << "DF(" << ccy << ",";
			dfDesc << "PayDate[i])";
		}
		rowDescVec[DF] = dfDesc.str();
		rowTypeVec[DF] = ARM_STRING;

		//FUNDING
		CC_Ostringstream fundingDesc;
		if (eventIdx == 0)
		{
			fundingDesc << "0.0";
		}
		else
		{
			fundingDesc << "(Libor[i-1]+Spread[i])*FundingBasis[i]*DF[i]*Notional[i]";
		}
		rowDescVec[Funding] = fundingDesc.str();
		rowTypeVec[Funding] = ARM_STRING;

		//RATIO ADV
		CC_Ostringstream ratioAdvDesc;
		if (eventIdx == 0)
		{
			ratioAdvDesc << "0.0";
		}
		else
		{
			ratioAdvDesc << "if(Index[i-1]<Strike[i]+";
			ratioAdvDesc << itsSnowRangeParams->Elt(0)/100.0; //barrier shift
			ratioAdvDesc << ",1,0)";
		}
		rowDescVec[RatioADV] = ratioAdvDesc.str();
		rowTypeVec[RatioADV] = ARM_STRING;

		//RATIO ARR
		CC_Ostringstream ratioArrDesc;
		if (eventIdx == 0)
		{
			ratioArrDesc << "0.0";
		}
		else
		{
			ratioArrDesc << "if(Index[i]<Strike[i]+";
			ratioArrDesc << itsSnowRangeParams->Elt(0)/100.0; //barrier shift
			ratioArrDesc << ",1,0)";
		}
		rowDescVec[RatioARR] = ratioArrDesc.str();
		rowTypeVec[RatioARR] = ARM_STRING;

		//NB PERIOD
		CC_Ostringstream nbPeriodDesc;
		if (eventIdx == 0)
		{
			nbPeriodDesc << "0.0";
		}
		else
		{
			double basis = CountYears(itsFundingDayCount, startDate, endDate);
			nbPeriodDesc << floor(basis/StringMatuToYearTerm(const_cast<char*>(itsCouponIndexTerm.c_str()))+0.5);
		}
		rowDescVec[NbPeriod] = nbPeriodDesc.str();
		rowTypeVec[NbPeriod] = ARM_DOUBLE;

		//FCORR
		CC_Ostringstream fCorrDesc;
		if (eventIdx == 0)
		{
			fCorrDesc << "0.0";
		}
		else
		{
			fCorrDesc << itsSnowRangeParams->Elt(3); //theta
			fCorrDesc << "*(RatioADV[i]+RatioARR[i])/NbPeriod[i]";
		}
		rowDescVec[FCorr] = fCorrDesc.str();
		rowTypeVec[FCorr] = ARM_STRING;

		//COUPON LEV
		CC_Ostringstream couponLevDesc;
		if (eventIdx == 0)
		{
			couponLevDesc << "0.0";
		}
		else if (eventIdx == 1)
		{
			couponLevDesc << "FCorr[i]";
		}
		else
		{
			couponLevDesc << "if(Ratchet[i]==1,CouponLev[i-1],1)*FCorr[i]";
		}
		rowDescVec[CouponLev] = couponLevDesc.str();
		rowTypeVec[CouponLev] = ARM_STRING;

		//COUPON
		CC_Ostringstream couponDesc;
		if (eventIdx == 0)
		{
			couponDesc << "0.0";
		}
		else
		{
			couponDesc << "-(" << const_cast< ARM_ReferenceValue& >(itsFixedRate).Interpolate(startDate);
			couponDesc << "+" << const_cast< ARM_ReferenceValue& >(itsLeverage).Interpolate(startDate);
			couponDesc << "*Index[i-1])*CouponLev[i]*DF[i]*CouponBasis[i]*Notional[i]";
		}
		rowDescVec[Coupon] = couponDesc.str();
		rowTypeVec[Coupon] = ARM_STRING;

		//STD COUPON
		CC_Ostringstream stdCouponDesc;
		if (eventIdx == 0)
		{
			stdCouponDesc << "0.0";
		}
		else
		{
			stdCouponDesc << "-(" << const_cast< ARM_ReferenceValue& >(itsFixedRate).Interpolate(startDate);
			stdCouponDesc << "+" << const_cast< ARM_ReferenceValue& >(itsLeverage).Interpolate(startDate);
			stdCouponDesc << "*Index[i-1])*FCorr[i]*DF[i]*CouponBasis[i]*Notional[i]";
		}
		rowDescVec[StdCoupon] = stdCouponDesc.str();
		rowTypeVec[StdCoupon] = ARM_STRING;

		//PLAIN COUPON
		CC_Ostringstream plainCouponDesc;
		if (eventIdx == 0)
		{
			plainCouponDesc << "0.0";
		}
		else
		{
			plainCouponDesc << "-(" << const_cast< ARM_ReferenceValue& >(itsFixedRate).Interpolate(startDate);
			plainCouponDesc << "+" << const_cast< ARM_ReferenceValue& >(itsLeverage).Interpolate(startDate);
			plainCouponDesc << "*Index[i-1])*DF[i]*CouponBasis[i]*Notional[i]";
		}
		rowDescVec[PlainCoupon] = plainCouponDesc.str();
		rowTypeVec[PlainCoupon] = ARM_STRING;

		//SWAP
		CC_Ostringstream swapDesc;
		if (eventIdx == 0)
		{
			swapDesc << "0.0";
		}
		else
		{
			swapDesc << "Funding[i]+Coupon[i]";
		}
		rowDescVec[Swap] = swapDesc.str();
		rowTypeVec[Swap] = ARM_STRING;

		//CALL OPTION
		CC_Ostringstream callOptionDesc;
		if (eventIdx == 0)
		{
			callOptionDesc << "0.0";
		}
		else if (eventIdx == nbEvent-1)
		{
			callOptionDesc << "Swap[i]";
		}
		else
		{
			callOptionDesc << "Exercise(Swap[i],-CF[i]*Notional[i],CallOption[i+1],1,Swap[i],Num[i],Swap[i]*Swap[i],Num[i]*Swap[i],Num[i]*Num[i],Swap[i]*Swap[i]*Swap[i],Num[i]*Swap[i]*Swap[i],Num[i]*Num[i]*Swap[i],Num[i]*Num[i]*Num[i])";
		}
		rowDescVec[CallOption] = callOptionDesc.str();
		rowTypeVec[CallOption] = ARM_STRING;

		//CALL SWAP
		CC_Ostringstream callSwapDesc;
		if (eventIdx == 1)
		{
			callSwapDesc << "CallOption[i]";
		}
		else
		{
			callSwapDesc << "0.0";
		}
		rowDescVec[CallSwap] = callSwapDesc.str();
		rowTypeVec[CallSwap] = ARM_STRING;

		return ARM_RowInfo(rowDescVec,rowTypeVec);
	}
	catch (Exception& x)
	{
	    x.DebugPrint();

		throw x;
	}
	catch (...)
	{
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in SnowRangeCalculator::MiddleRows" );
	}
}

/*****************************************************************************
/ Class  : ARM_SnowRangeCalculator
/ Routine: CreateAndSetModel
/ Returns: void
/ Action : create Markov Fonctional or Smiled FRM
*****************************************************************************/
void ARM_SnowRangeCalculator::CreateAndSetModel()
{
	/// Create MC Method

	// random generator
	ARM_RandomGenerator* randGen1=NULL;
	ARM_RandomGenerator* randGen2=NULL;
	ARM_RandomGenerator* randGen=NULL;

	ARM_RandomGeneratorPtrVector randGenVector;

	ARM_RandomGeneratorPtr pRandGen1, pRandGen2;
		
	randGen1 = ARM_RandGenFactory.Instance()->CreateRandGen(
			(ARM_RandGenFactoryImp::BaseGenType) ARM_ArgConv_BaseGenAlgoType.GetNumber(itsGeneratorType));
	pRandGen1 = ARM_RandomGeneratorPtr(randGen1);

	randGen2 = ARM_RandGenFactory.Instance()->CreateRandGen(
			ARM_RandGenFactoryImp::UnknownBaseGenAlgorithm,
			(ARM_RandGenFactoryImp::TransformAlgo) ARM_ArgConv_TransformAlgoType.GetNumber(itsInversionMethod),
			pRandGen1);

	if (itsAntithetic)
	{
		pRandGen2 = ARM_RandomGeneratorPtr(randGen2);
		randGen = ARM_RandGenFactory.Instance()->CreateRandGen(
			ARM_RandGenFactoryImp::UnknownBaseGenAlgorithm,
			(ARM_RandGenFactoryImp::TransformAlgo) ARM_ArgConv_TransformAlgoType.GetNumber("AntitheticOne"),
			pRandGen2);
	}
	else
		randGen = randGen2;

	randGenVector.resize(1);
	randGenVector[0] = ARM_RandomGeneratorPtr(randGen);

	// scheduler
	ARM_SchedulerBase* scheduler = NULL;
	ARM_SamplerBase* sampler = NULL;
	int multiDim = 2;
	std::vector<double> samplerDatas(0);

	int schedulerType = ARM_SchedulerBase::TimeStepPerYear;
	std::vector<double> schedulerDatas(1);
	schedulerDatas[0] = 1;

	int nbMCSteps;
	scheduler = ARM_SchedulerFactory.Instance()->CreateScheduler(schedulerType,
																 schedulerDatas,
																 nbMCSteps);

	sampler = ARM_SamplerFactory.Instance()->CreateSampler( multiDim,
															itsSamplerType,
															samplerDatas,
															scheduler );

	ARM_AMCLongstaffSchwartz * AMCLongstaffSchwartz = new ARM_AMCLongstaffSchwartz( itsNbSteps );
	
	// American Monte Carlo method
	// it will be set to model after calibration !
	itsMCMethod = ARM_NumMethodPtr(new ARM_AMCMethod(itsNbSteps,
													 randGenVector,
													 sampler,
													 AMCLongstaffSchwartz));

	ARM_ZeroCurve* zc = dynamic_cast< ARM_BSSmiledModel* >( GetMktDataManager()->GetData(GetKeys()[CfModelKey]) )->GetZeroCurve();

	// Create Model
	if (itsModelName == "HK")
	{
		// HK parameters
		ARM_CurveModelParam* meanReversion = dynamic_cast< ARM_CurveModelParam* >( GetMktDataManager()->GetData(GetKeys()[MrsKey]) );

		int size = (*itsProductDateStrip).GetResetDates()->size();
		std::vector<double> breakPointTimes(size);
		std::vector<double> values(size);

		ARM_Date AsOf = zc->GetAsOfDate();

		for (int i=0; i<size; i++)
		{
			breakPointTimes[i] = (*itsProductDateStrip).GetResetDates()->Elt(i) - AsOf.GetJulian();
			values[i] = itsModelParams->Elt(0) + i*itsModelParams->Elt(1);
		}

		ARM_CurveModelParam calibParamHKVol = ARM_CurveModelParam(ARM_ModelParamType::Volatility, 
																	 &values, 
																	 &breakPointTimes);

		vector<ARM_ModelParam*> paramsVector;
		paramsVector.push_back((ARM_ModelParam*)meanReversion);
		paramsVector.push_back((ARM_ModelParam*)&calibParamHKVol);

		// Create MF Model
		ARM_MarkovFunctional* mod = new ARM_MarkovFunctional( CreateClonedPtr( zc ) );
		
		ARM_ModelParamVector modelParams(paramsVector);
		mod->SetModelParams(ARM_ModelParamsMF(modelParams));

		ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::TerminalZc ) );
		mod->SetNumeraire( numeraire );

		SetPricingModel(ARM_PricingModelPtr(static_cast< ARM_PricingModel* >(mod)));
	}
	else if (itsModelName == "SBGM")
	{
		std::vector<double> breakPointTimes(1);
		std::vector<double> values(1);
		std::vector<double> lowerBound(1);
		std::vector<double> upperBound(1);
		
		breakPointTimes[0] = 0;
		values[0] = itsCalibParams->Elt(3) ? 0.2 : itsModelParams->Elt(0)/100.0;
		lowerBound[0] = 0;
		upperBound[0] = 1;

		ARM_CurveModelParam betaCorrelParam = ARM_CurveModelParam ( ARM_ModelParamType::BetaCorrelation, 
																	&values, 
																	&breakPointTimes,
																	"",
																	"STEPUPRIGHT",
																	&lowerBound,
																	&upperBound );

		// Create Smiled FRM Model
		int timeStepsNb = itsCalibParams->Elt(0);
		int gridSize = itsCalibParams->Elt(1);
		int stdDevNb = itsCalibParams->Elt(2);
		bool skipPDE = false;
		bool allowInterpol = false;
		string swaptionApprox("LOCAL+");
		
		ARM_SmiledFRM* mod = new ARM_SmiledFRM( CreateClonedPtr( zc ),
												NULL,
												timeStepsNb,
												gridSize,
												stdDevNb,
												skipPDE,
												allowInterpol,
												(ARM_ModelParamsSmiled::CalibProxy) ARM_ArgConv_MMCalibProxy.GetNumber(swaptionApprox) );

		ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::TerminalZc ) );
		mod->SetNumeraire( numeraire );

		vector<ARM_ModelParam*> paramsVector;
		paramsVector.push_back((ARM_ModelParam*)&betaCorrelParam);
		ARM_CurveModelParam* hump = dynamic_cast< ARM_CurveModelParam* >( GetMktDataManager()->GetData(GetKeys()[HumpKey]) );
		paramsVector.push_back((ARM_ModelParam*)hump);

		int factorsNb = itsModelParams->Elt(1);
		mod->SetModelParams(ARM_ModelParamsSmiled(paramsVector, factorsNb));

		SetPricingModel(ARM_PricingModelPtr(static_cast< ARM_PricingModel* >(mod)));
	}
}

/*****************************************************************************
/ Class  : ARM_GlobaCapCalculator
/ Routine: CreateAndSetCalibration
/ Returns: void
/ Action : create calibMethod
*****************************************************************************/
void ARM_SnowRangeCalculator::CreateAndSetCalibration()
{
	try
	{
		ARM_StdPortfolio* swoptPf = NULL;
		ARM_CalibMethod* calibMethod = NULL;

		int size = (*itsCalibDateStrip).GetResetDates()->size();

		// Create a Numerical Calib Method without portfolio
		// For HK : with a portfolio of swaptions if we need to calibrate swaptions, without portfolio otherwise.
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

		// build another date strip for calibration schedule
		// when fixing is ARREARS, we should deal with FORWARD START and FORWARD END dates
		// instead of start and end dates !
		ARM_DateStrip* dateStrip = (ARM_DateStrip*)itsCalibDateStrip->Clone();
		if (itsResetTiming == K_ARREARS)
		{
			dateStrip->SetFlowStartDates(dateStrip->GetFwdRateStartDates());
			dateStrip->SetFlowEndDates(dateStrip->GetFwdRateEndDates());
		}

		for (int i = 0; i < size; ++i)
		{
			ARM_VanillaSecurityDensity* vanillaSecDensity	= NULL;
			double matu = ((*(dateStrip->GetResetDates()))[i] - AsOf.GetJulian())/365.0;
			double tenor = ((*(dateStrip->GetFlowEndDates()))[i] - (*(dateStrip->GetFlowStartDates()))[i])/365.0;
			
			double alpha = alphaCurve->ComputeVolatility(matu, tenor)/100.0;
			double rho = rhoCurve->ComputeVolatility(matu, tenor);
			double nu = nuCurve->ComputeVolatility(matu, tenor);
			double beta = betaCurve->ComputeVolatility(matu, tenor);

			ARM_SABRDensityFunctor* method = new ARM_SABRDensityFunctor(alpha, beta, rho, nu, capModel->GetSABRFlag());
			
			int frequency = K_ANNUAL;
			if ((itsCcy.GetCcyName() == "USD") && (StringMatuToYearTerm(const_cast<char*>(itsCouponIndexTerm.c_str())) > 1.0))
				frequency = K_SEMIANNUAL;

			vanillaSecDensity = new ARM_VanillaSecurityDensity( (*(dateStrip->GetResetDates()))[i], 
																(*(dateStrip->GetFlowStartDates()))[i],
																(*(dateStrip->GetFlowEndDates()))[i],
																ARM_DensityFunctorPtr(method),
																frequency ); 

			vanillaSecDensities[i] = ARM_VanillaSecDensityPtr( static_cast<ARM_VanillaSecurityDensity*> (vanillaSecDensity->Clone()) );
		}

		if (itsModelName == "HK" && itsCalibParams->Elt(3)) //if calibrate swaption, create a portfolio of swaptions
		{
			ARM_BSSmiledModel* swoptModel = dynamic_cast< ARM_BSSmiledModel* >( GetMktDataManager()->GetData(GetKeys()[OswModelKey]) );
			int indexType = FromFrequencyToXiborType(itsPayFreq, zc->GetCurrencyUnit()->GetCcyName());
			
			double maturityDate = (*(dateStrip->GetFwdRateEndDates()))[size-1];
		
			swoptPf = new ARM_StdPortfolio();
			
			for (i = 1; i < size-1; i++) // HK case : n-2 swaptions
			{
				double startJulian = (*(dateStrip->GetFwdRateStartDates()))[i];

				ARM_Swaption* swaption = new ARM_Swaption(ARM_Date(startJulian), 
														  ARM_Date(maturityDate), 
														  K_RCV, 
														  K_EUROPEAN, 
														  K_MARKET_RATE,
														  ARM_Date(startJulian).PreviousBusinessDay(-itsResetGap, const_cast<char*>(itsResetCal.c_str())),
														  (ARM_INDEX_TYPE)indexType,
														  0.0,
														  -1000000.0,
														  K_DEF_FREQ,
														  K_DEF_FREQ,
														  &itsCcy);

				swaption->SetModelVariable(NULL);
				swaption->SetModel(swoptModel);

				double mktPrice = swaption->ComputePrice();
				double curWeight = ( itsCashFlow.GetDiscreteValues()->Elt(i) > 10 ) ? 0.0 : 1.0;       

				swoptPf->AddInstrument(swaption, mktPrice, curWeight, 0);
			}
		}
		
		ARM_ModelParamVector calibParams = ARM_ModelParamVector(0);
		ARM_ModelParams* Modelparams = &*GetPricingModel()->GetModelParams();
		if (Modelparams)
		{
			ARM_ModelParamVector modelParamsVect = Modelparams->GetModelParams();
			calibParams = ARM_ModelParamVector(1, modelParamsVect[ARM_ModelParamType::Volatility]);
		}

		calibMethod = new ARM_CalibMethod ( ARM_StdPortfolioPtr(swoptPf),
											calibParams,
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
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in SnowRangeCalculator::CreateAndSetCalibration" );
	}
}

/*****************************************************************************
/ Class  : ARM_GlobaCapCalculator
/ Routine: CreateCalibMethodForSBGM
/ Returns: void
/ Action : Create calibMethod
/          This fonction is called only if we price with SBGM
/          and calibrate swaptions
*****************************************************************************/
void ARM_SnowRangeCalculator::CreateCalibMethodForSBGM()
{
	try
	{
		ARM_StdPortfolio* swoptPf = new ARM_StdPortfolio();
		ARM_CalibMethod* calibMethod = NULL;
		
		int size = (*itsCalibDateStrip).GetResetDates()->size();

		ARM_BSSmiledModel* swoptModel = dynamic_cast< ARM_BSSmiledModel* >( GetMktDataManager()->GetData(GetKeys()[OswModelKey]) );
		ARM_ZeroCurve* zc = swoptModel->GetZeroCurve();
		int indexType = FromFrequencyToXiborType(itsPayFreq, zc->GetCurrencyUnit()->GetCcyName());
		
		double maturityDate = (*(itsCalibDateStrip->GetFwdRateEndDates()))[size-1];
	
		// create a portfolio of swaptions
		for (int i=1; i<size; i++) // SBGM case : n-1 swaptions
		{
			double startJulian = (*(itsCalibDateStrip->GetFwdRateStartDates()))[i];

			ARM_Swaption* swaption = new ARM_Swaption(ARM_Date(startJulian), 
													  ARM_Date(maturityDate), 
													  K_RCV, 
													  K_EUROPEAN, 
													  K_MARKET_RATE,
													  ARM_Date(startJulian).PreviousBusinessDay(-itsResetGap, const_cast<char*>(itsResetCal.c_str())),
													  (ARM_INDEX_TYPE)indexType,
													  0.0,
													  -1000000.0,
													  K_DEF_FREQ,
													  K_DEF_FREQ,
													  &itsCcy);

			swaption->SetModelVariable(NULL);
			swaption->SetModel(swoptModel);

			double mktPrice = swaption->ComputePrice();
			double curWeight = ( itsCashFlow.Interpolate(startJulian) > 10 ) ? 0.0 : 1.0;       

			swoptPf->AddInstrument(swaption, mktPrice, curWeight, 0);
		}

		ARM_ModelParams* Modelparams = &*GetPricingModel()->GetModelParams();
		ARM_ModelParamVector modelParamsVect = Modelparams->GetModelParams();
		ARM_ModelParamVector calibParams(1, modelParamsVect[ARM_ModelParamType::BetaCorrelation]);

		///Create Calib Method
		calibMethod = new ARM_CalibMethod(ARM_StdPortfolioPtr(swoptPf),
										  calibParams,
										  ARM_CalibMethodType::Bootstrap1D, 
										  100, 
										  ARM_CalibrationTarget::PriceTarget,
										  NULL, 
										  NULL, 
										  false);

		SetCalibMethod(ARM_CalibMethodPtr(calibMethod));
	}
	catch (Exception& x)
	{
	    x.DebugPrint();

		throw x;
	}
	catch (...)
	{
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in SnowRangeCalculator::CreateCalibMethodForSBGM" );
	}
}

void ARM_SnowRangeCalculator::UpdateModel()
{
	CreateAndSetModel();
}

void ARM_SnowRangeCalculator::UpdateCalibration(bool isUpdateStrike)
{
	CreateAndSetCalibration();
}

/*****************************************************************************
/ Class  : ARM_GlobaCapCalculator
/ Routine: Calibrate
/ Returns: void
/ Action : Calibrate the model, set MC Method afterwards
*****************************************************************************/
void ARM_SnowRangeCalculator::Calibrate()
{
	try
	{
		if (itsModelName == "HK")
		{
			// Create PDE Method
			ARM_PDEMethod* PDEMethod= NULL;

			// Gets the right NumScheme from its name
			int numSchemeType = ARM_ArgConv_PDENumSchemeType.GetNumber(string("CN1F"));
			ARM_PDENumericalScheme* numScheme = ARM_PDENumericalScheme::getNumericalSchemeInstanceById(numSchemeType);

			// Builds PDE NumMethod With this numscheme
			PDEMethod = new ARM_PDEMethod(numScheme, itsCalibParams->Elt(0), itsCalibParams->Elt(1), 2.*itsCalibParams->Elt(2));
			(&*GetPricingModel())->SetNumMethod(ARM_NumMethodPtr(PDEMethod));
		}

		GetCalibMethod()->Calibrate(&*GetPricingModel());

		// if we price with a SBGM model and want to calibrate swaptions, 
		// we'll need to RECALIBRATE using a calib method containing a portfolio of swaptions
		if (itsModelName == "SBGM" && itsCalibParams->Elt(3))
		{
			CreateCalibMethodForSBGM();
			GetCalibMethod()->Calibrate(&*GetPricingModel());
		}
		
		// Set MC Method for pricing
		(&*GetPricingModel())->SetNumMethod(itsMCMethod);
	}
	catch (Exception& x)
	{
	    x.DebugPrint();

		throw x;
	}
	catch (...)
	{
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in SnowRangeCalculator::Calibrate" );
	}
}

	/// Initialisation to 0 of all columns of the deal description that could be priced
void ARM_SnowRangeCalculator::InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const
{
    string zeroValue("0");

	rowDescVec[CallSwap] = zeroValue;
	rowTypeVec[CallSwap] = ARM_DOUBLE;

	rowDescVec[CallOption] = zeroValue;
	rowTypeVec[CallOption] = ARM_DOUBLE;

	rowDescVec[Swap] = zeroValue;
	rowTypeVec[Swap] = ARM_DOUBLE;

	rowDescVec[PlainCoupon] = zeroValue;
	rowTypeVec[PlainCoupon] = ARM_DOUBLE;

	rowDescVec[StdCoupon] = zeroValue;
	rowTypeVec[StdCoupon] = ARM_DOUBLE;

	rowDescVec[Coupon] = zeroValue;
	rowTypeVec[Coupon] = ARM_DOUBLE;

	rowDescVec[Funding] = zeroValue;
	rowTypeVec[Funding] = ARM_DOUBLE;

	rowDescVec[Libor] = zeroValue;
	rowTypeVec[Libor] = ARM_DOUBLE;

	rowDescVec[Index] = zeroValue;
	rowTypeVec[Index] = ARM_DOUBLE;

	rowDescVec[DF] = zeroValue;
	rowTypeVec[DF] = ARM_DOUBLE;

	rowDescVec[RatioADV] = zeroValue;
	rowTypeVec[RatioADV] = ARM_DOUBLE;

	rowDescVec[RatioARR] = zeroValue;
	rowTypeVec[RatioARR] = ARM_DOUBLE;

	rowDescVec[FCorr] = zeroValue;
	rowTypeVec[FCorr] = ARM_DOUBLE;
}

/*****************************************************************************
/ Class  : ARM_SnowRangeCalculator
/ Routine: PricedColumnNames
/ Returns: ARM_StringVector
/ Action : create the priced column names of the deal description
******************************************************************************/
ARM_StringVector ARM_SnowRangeCalculator::PricedColumnNames() const
{
	ARM_StringVector pricedColumns;

	// effective nb of columns to price :
	int colSize = MIN(NbProductsToPrice, itsProductsToPrice.size());

	for (int i=0; i<colSize; i++)
	{
		if (itsProductsToPrice[i])
			pricedColumns.push_back(SnowRangeColNamesTable[SnowRangeProductToPriceColumns[i]]);
	}

	return pricedColumns;
}


void ARM_SnowRangeCalculator::ComputePricingData() const
{
	if (!itsHasBeenPriced)
		const_cast<ARM_SnowRangeCalculator*>(this)->PriceAndTimeIt();

	GetPricingData()[ "CallSwap" ]			= itsCallSwapPrice;
	GetPricingData()[ "CallOption" ]		= itsCallOptionPrice;
	GetPricingData()[ "Swap" ]				= itsSwapPrice;
	GetPricingData()[ "PlainCoupon" ]		= itsPlainCouponPrice;	
	GetPricingData()[ "StdCoupon" ]			= itsStdCouponPrice;	
	GetPricingData()[ "Coupon" ]			= itsCouponPrice;
	GetPricingData()[ "Funding" ]			= itsFundingPrice;
	GetPricingData()[ "Libor" ]				= itsLiborPrice;
	GetPricingData()[ "Index" ]				= itsIndexPrice;
	GetPricingData()[ "DF" ]				= itsDFPrice;
	GetPricingData()[ "RatioADV" ]			= itsRatioADVPrice;
	GetPricingData()[ "RatioARR" ]			= itsRatioARRPrice;
	GetPricingData()[ "FCorr" ]				= itsFCorrPrice;
	GetPricingData()[ "CouponLev" ]			= itsCouponLevPrice;
}

/*****************************************************************************
/ Class  : ARM_SnowRangeCalculator
/ Routine: Price
/ Returns: the call swap value
/ Action : price the call
*****************************************************************************/
double ARM_SnowRangeCalculator::Price()
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
				price	= genPricer->GetPricerInfo()->GetContents( SnowRangeColNamesTable[ SnowRangeProductToPriceColumns[i] ] ).GetData("Price").GetDouble();
						
				if (i == CallSwapPrice)
				{
					itsCallSwapPrice = price;
				}
				else if (i == CallOptionPrice)
				{
					itsCallOptionPrice = price;
				}
				else if (i == SwapPrice)
				{
					itsSwapPrice = price;
				}
				else if (i == PlainCouponPrice)
				{
					itsPlainCouponPrice = price;
				}
				else if (i == StdCouponPrice)
				{
					itsStdCouponPrice = price;
				}
				else if (i == CouponPrice)
				{
					itsCouponPrice = price;
				}
				else if (i == FundingPrice)
				{
					itsFundingPrice = price;
				}
				else if (i == LiborPrice)
				{
					itsLiborPrice = price;
				}
				else if (i == IndexPrice)
				{
					itsIndexPrice = price;
				}
				else if (i == DFPrice)
				{
					itsDFPrice = price;
				}
				else if (i == RatioADVPrice)
				{
					itsRatioADVPrice = price;
				}
				else if (i == RatioARRPrice)
				{
					itsRatioARRPrice = price;
				}
				else if (i == FCorrPrice)
				{
					itsFCorrPrice = price;
				}
				else if (i == CouponLevPrice)
				{
					itsCouponLevPrice = price;
				}
			}
		}
		itsHasBeenPriced = true;

		price = itsCallSwapPrice * itsPayRec;
		return price;
	}
	catch (Exception& x)
	{
	    x.DebugPrint();

		throw x;
	}
	catch (...)
	{
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in SnowRangeCalculator::Price" );
	}
}

/*****************************************************************************
/ Class   : ARM_SnowRangeCalculator
/ Routines: GeneralDataToString
/ Returns :
/ Action  : Construct information string for
/           general product data
******************************************************************************/
string ARM_SnowRangeCalculator::GeneralDataToString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream generalGCData;

	//General Data
	generalGCData << "\nGENERAL DATA\n"	<< endl;
	generalGCData << "Currency\t\t: "	<<  itsCcy.GetCcyName() << endl;
	generalGCData << "StartDate\t\t: "  <<  itsStartDate.toString() << endl;
    generalGCData << "EndDate\t\t\t: "	<<  itsEndDate.toString() << endl;
    generalGCData << "Pay/Rec\t\t\t: "	<<  ARM_ParamView::GetMappingName(S_RECEIVE_PAY, itsPayRec) << endl;
	generalGCData << "\n";
	
	//Funding Data
	generalGCData << "FUNDING DATA\n"		<< endl;
	generalGCData << "Index Term\t\t: "		<< itsFundingIndexTerm << endl;
	generalGCData << "Day Count\t\t: "		<< ARM_ArgConvReverse_DayCount.GetString(itsFundingDayCount) << endl;

	//Coupon Data
	generalGCData << "COUPON DATA\n"			<< endl;
	generalGCData << "Index Term\t\t: "			<< itsCouponIndexTerm << endl;
	generalGCData << "Day Count\t\t: "			<< ARM_ArgConvReverse_DayCount.GetString(itsCouponDayCount) << endl;
	generalGCData << "Reset Frequency\t\t: "	<< ARM_ArgConvReverse_StdFrequency.GetString(itsResetFreq) << endl;
	generalGCData << "Pay Frequency\t\t: "		<< ARM_ArgConvReverse_StdFrequency.GetString(itsPayFreq) << endl;
	generalGCData << "Reset Timing\t\t: "		<< ARM_ArgConvReverse_Timing.GetString(itsResetTiming) << endl;
	generalGCData << "Pay Timing\t\t: "			<< ARM_ArgConvReverse_Timing.GetString(itsPayTiming) << endl;
	generalGCData << "Reset Calendar\t\t: "		<< itsResetCal << endl;
	generalGCData << "Pay Calendar\t\t: "		<< itsPayCal << endl;
	generalGCData << "\n";

	generalGCData << endl;

    return generalGCData.str();
}

/*****************************************************************************
/ Class   : ARM_SnowRangeCalculator
/ Routines: DealDesDataToString
/ Returns :
/ Action  : Construct information string for
/			deal description data
******************************************************************************/
string ARM_SnowRangeCalculator::DealDesDataToString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream dealDesGCData;

	double asOf = GetMktDataManager()->GetAsOfDate().GetJulian();
	//Notional Profile
	dealDesGCData << "=====================Notional Profile===================" << endl;
	ARM_GP_CurvePtr notional(RefValueToCurve(itsNotional, asOf));
	dealDesGCData << notional->toString();

	//Spread Profile
	dealDesGCData << "\n=====================Spread Profile===================" << endl;
	ARM_GP_CurvePtr spread(RefValueToCurve(itsSpread, asOf));
	dealDesGCData << spread->toString();

	//Strike Profile
	dealDesGCData << "\n======================Strike Profile===================" << endl;
	ARM_GP_CurvePtr strike(RefValueToCurve(itsStrike, asOf));
	dealDesGCData << strike->toString();

	//Ratchet Profile
	dealDesGCData << "\n=====================Ratchet Profile==================" << endl;
	ARM_GP_CurvePtr ratchet(RefValueToCurve(itsRatchet, asOf));
	dealDesGCData << ratchet->toString();

	//Cash Flow Profile
	dealDesGCData << "\n=====================CashFlow Profile=================" << endl;
	ARM_GP_CurvePtr cashFlow(RefValueToCurve(itsCashFlow, asOf));
	dealDesGCData << cashFlow->toString();

	//Fixed Rate Profile
	dealDesGCData << "\n===================Fixed Rate Profile=================" << endl;
	ARM_GP_CurvePtr fixedRate(RefValueToCurve(itsFixedRate, asOf));
	dealDesGCData << fixedRate->toString();

	//Leverage Profile
	dealDesGCData << "\n===================Leverage Profile=================" << endl;
	ARM_GP_CurvePtr leverage(RefValueToCurve(itsLeverage, asOf));
	dealDesGCData << leverage->toString();

    return dealDesGCData.str();
}

/*****************************************************************************
/ Class  : ARM_SnowRangeCalculator
/ Routine: View
/ Returns: 
/ Action : .
******************************************************************************/
void ARM_SnowRangeCalculator::View(char* id, FILE* ficOut) const
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
	fprintf(fOut,"\n===============> SNOW RANGE CALCULATOR <================");
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
