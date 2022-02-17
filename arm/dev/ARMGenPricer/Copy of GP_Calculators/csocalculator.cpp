/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file csocalculator.cpp
 *
 *  \brief file for the calculator for Callable SpreadOption
 *
 *	\author  JP Riaudel
 *	\version 1.0
 *	\date May 2005
 */

#include "gpbase/removeidentifiedwarning.h"
#include "gpcalculators/csocalculator.h"
#include "gpcalculators/localcsocalculator.h"

/// gpbase
#include "gpbase/env.h"
#include "gpbase/singleton.h"
#include "gpbase/autocleaner.h"
#include "gpbase/datestripcombiner.h"
#include "gpbase/gpmatrixlinalg.h"

/// gpinfra
#include "gpinfra/pricingmodel.h"
#include "gpinfra/pricingadviser.h"
#include "gpinfra/gensecurity.h"
#include "gpinfra/genpricer.h"
#include "gpinfra/dealdescription.h"
#include "gpinfra/numeraire.h"
#include "gpinfra/numerairefactory.h"
#include "gpinfra/correlmatparam.h"

/// gpmodels
#include "gpmodels/modelparamssfrm.h"
#include "gpmodels/modelparamssfrmfactory.h"
#include "gpmodels/sfrm.h"

/// gpcalib
#include "gpcalib/calibmethod.h"
#include "gpcalib/vanillaspreadoption.h"
#include "gpcalib/kerneltogp.h"
#include "gpcalib/stripper.h"

/// gpnummethods
#include "gpnummethods/treebase.h"
#include "gpnummethods/treefactory.h"
#include "gpnummethods/scheduler.h"
#include "gpnummethods/sampler.h"
#include "gpnummethods/truncator.h"
#include "gpnummethods/reconnector.h"
#include "gpnummethods/smoother.h"

#include <crv/zerocurv.h>
#include <crv/volflat.h>
#include <util/fromto.h>
#include <inst/swaption.h>
#include <inst/spreadoption.h>
#include <inst/forex.h>
#include <crv/correlmanager.h>


CC_BEGIN_NAMESPACE( ARM )


// SFRM default factors number
const int SFRM_NB_FACTORS			= 2;
const int SFRM_VOL_TYPE				= K_ROW;

const int TREE_NBSTEPS_PER_YEAR	=5;
const double STD_DEV_RATIO		=5.0;
const double MIN_STD_DEV		=0.001;

const double DEFAULT_PRICE=1.0e+100;
const double DEFAULT_PRECISION=0.000001;
const double DEFAULT_WEIGHT=1.0;




const string ARM_CSOCalculator::CSOColNamesTable [] =
{
    "ResetDate",
    "StartDate",
    "PayDate",
    "EndDateCMS1",
    "EndDateCMS2",
    "IT",
    "CMS1",
    "CMS2",
    "CpnFix",
	"Strike",
	"Notional",
	"Leverage",
    "SpreadOptionLet",
	"SpreadOption",
	"StrikeFlow",
	"PVStrikeFlow",
    "FundStartDate",
    "FundEndDate",
	"FundSpread",
	"FundFragmentNPV",
	"FundNPV",
    "Flow",
    "TotalFlow",
    "Fees",
	"BermudaProfile",
	"CSO"/*,
	"ExerciseCondition",
	"NumeraireDate",
	"ExerciseConditionOfIndex",
	"ProbaOfExercise",
	"Funding",
	"Strike1",
	"Strike2",
	"CapCMS1",
	"CapCMS2"*/
};

/// Reference schedules for CSO date structure
const unsigned int EXERCISE_SCHED =0;
const unsigned int NB_CSO_SCHED =1;


/// SFRM sigma range [1bp,10000bp] with a 500bp default value
const double SIGMA_LOWER_BOUND      = 0.0001;
const double SIGMA_UPPER_BOUND      = 1.0;
const double SIGMA_DEFAULT_VALUE    = 0.20;
const int SIGMA_CSO_ARM_MAX_ITER = 5;
const double SIGMA_DEFAULT_PRECISION= 0.000001;

/// SFRM beta range [0%,130%] with a 100% default value
const double BETA_DEFAULT_VALUE    = 1;
const double BETA_LOWER_BOUND	   = 0.01;
const double BETA_UPPER_BOUND	   = 1.3;
const int BETA_CSO_ARM_MAX_ITER = 5;
const double BETA_DEFAULT_PRECISION= 0.000001;

/// SFRM MRS range [-15%,50%] with a 2% default value
const double MRS_LOWER_BOUND        = -0.15;
const double MRS_UPPER_BOUND        = 0.15;
const double MRS_DEFAULT_VALUE      = 0.0;
const int MRS_CSO_ARM_MAX_ITER	= 5;
const double MRS_DEFAULT_PRECISION= 0.00001;

const double CORREL_CSO_ARM_MAX_ITER	= 5;
const double CORREL_DEFAULT_PRECISION= 0.00001;
const double CORREL_DEFAULT_PONDERATION = 0.1;
const double CORREL_DEFAULT_INITVALUE = 0.7;
const double CORREL_LOWER_BOUND      = -5.0;
const double CORREL_UPPER_BOUND      = 5.0;
const unsigned int CORREL_FREQ       = 1;

ARM_CSOCalculator::ARM_CSOCalculator(const ARM_Date& startDate,
									 const ARM_Date& endDate,
									 int CMSLong,
									 int CMSShort,
									 int cpnDayCount,
									 int cpnFreq,
									 int cpnResetTiming,
									 const ARM_Curve& nominal,
									 const ARM_Curve& fixCoupon,
									 const ARM_Curve& leverageLong,
									 const ARM_Curve& leverageShort,
									 const ARM_Curve& cpnMin,
									 const ARM_Curve& cpnMax,
									 const ARM_Curve& strike,
									 int fundFreq,
									 int fundDaycount,
									 const ARM_Curve& fundSpread,
									 int exerciseFreq,
									 int noticeGap,
									 int payRec,
									 const ARM_Curve& fees,
									 CalibrationMode sigmaCalibFlag,
									 CalibrationMode NDFlag,
									 bool mrsCalibFlag,
									 bool thetaCalibFlag,
									 std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice,
									 vector<double> ModelDatas,
									 const ARM_MarketData_ManagerRep& mktDataManager,
									 const ARM_StringVector& mdmKeys)
: ARM_GenCSOCalculator(startDate,
					   endDate,
					   CMSLong,
					   CMSShort,
					   cpnDayCount,
					   cpnFreq,
					   cpnResetTiming,
					   nominal,
					   fixCoupon,
					   leverageLong,
					   leverageShort,
					   cpnMin,
					   cpnMax,
					   strike,
					   fundFreq,
					   fundDaycount,
					   fundSpread,
					   ARM_FlatCurve(1.0),
					   exerciseFreq,
					   noticeGap,
					   payRec,
					   fees,
					   productsToPrice,
					   ModelDatas,
					   mktDataManager,
					   mdmKeys),
	itsSigmaCalib(sigmaCalibFlag),
	itsCalibND(NDFlag),
	itsMrsCalib(mrsCalibFlag),
	itsThetaCalib(thetaCalibFlag),
	itsModelDatas(ModelDatas),
	itsUpdateInitCalib(true)

{
	/// Check input datas
	CheckDataAndTimeIt();

	if (itsCalibND == TWODCalib)
		itsMrsCalib = false;

	/// Set the coupon & funding model name alias i.e. keys to access to
	/// the right model for each GP keyword (1st argument)
//    SetModelKeys();

	/// Create the Generic Security paid in coupon currency
	CreateAndSetDealDescriptionAndTimeIt();

	ARM_GP_Vector* resetDates = itsCalibDateStrip->GetResetDates();
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	int resetSize=resetDates->size();
	int i=0;
	while( (i < resetSize) && (const_cast< ARM_Curve& >(GetFixCpn()).Interpolate((*resetDates)[i]-asOfDate) > 0.0) )
		++i;
	SetNbFixFlows(i);

	/// Create the SFRM pricing model with its default parameters
	CreateAndSetModelAndTimeIt();

	/// Create the calibration set for volatility bootstapping and
	/// strike spread adjustments
	CreateAndSetCalibrationAndTimeIt();
}



ARM_CSOCalculator::ARM_CSOCalculator( const ARM_Date& asOfDate,
									  const ARM_Date& startDate,
									  const ARM_Date& endDate,
									  int CMSLong,
									  int CMSShort,
									  int cpnDayCount,
									  int cpnFreq,
									  int cpnResetTiming,
									  const ARM_Curve& nominal,
									  const ARM_Curve& fixCoupon,
									  const ARM_Curve& leverageLong,
									  const ARM_Curve& leverageShort,
									  const ARM_Curve& cpnMin,
									  const ARM_Curve& cpnMax,
									  const ARM_Curve& strike,
									  int fundFreq,
									  int fundDayCount,
									  const ARM_Curve& fundSpread,
									  int exerciseFreq,
									  int noticeGap,
									  int payRec,
									  const ARM_Curve& fees)
: ARM_GenCSOCalculator(asOfDate,
					   startDate,
                       startDate, // Would be fixEndDate
					   endDate,
					   CMSLong,
					   CMSShort,
					   cpnDayCount,
					   cpnFreq,
					   cpnResetTiming,
					   nominal,
					   fixCoupon,
					   leverageLong,
					   leverageShort,
					   cpnMin,
					   cpnMax,
					   strike,
					   fundFreq,
					   fundDayCount,
					   fundSpread,
					   ARM_FlatCurve(1.0),
					   exerciseFreq,
					   noticeGap,
					   payRec,
					   fees)
{
}

void ARM_CSOCalculator::InitCSOFromSummit( const ARM_MarketData_ManagerRep* mktDataManager,
									  const ARM_StringVector& mdmKeys,
									  CalibrationMode sigmaCalibFlag,
									  CalibrationMode NDFlag,
									  bool mrsCalibFlag,
									  bool thetaCalibFlag,
									  std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice,
									  vector<double> ModelDatas)
{
	InitializeMktDataManagerOnly(*mktDataManager);

	itsSigmaCalib = sigmaCalibFlag;
	itsCalibND = NDFlag;
	itsMrsCalib = mrsCalibFlag;
	itsThetaCalib = thetaCalibFlag;
	itsProductsToPrice = productsToPrice;
	itsModelDatas = ModelDatas;

    /// Set keys for MDM datas access
	///  We suppose that OSW smiled model  and OSW Model are same
    if(mdmKeys.size() < NbKeys)
    {
        /// To be compatible with non basis version
        ARM_StringVector newMdMKeys(mdmKeys);
        newMdMKeys.resize(NbKeys);
        newMdMKeys[FundingKey]		  = newMdMKeys[YcKey];
        newMdMKeys[BasisKey]		  = newMdMKeys[YcKey];
        newMdMKeys[ForexKey]		  = UNKNOWN_KEY_NAME;
        SetKeys(newMdMKeys);
    }
    else
        SetKeys(mdmKeys);

	/// Set the coupon/payment currency (inherited from ARM_Security)
    ARM_Currency* cpnCcy = static_cast< ARM_ZeroCurve* >(GetMktDataManager()->GetData(GetKeys()[YcKey]))->GetCurrencyUnit();
    SetCurrencyUnit(cpnCcy);

    /// Set funding & basis currencies
    ARM_Currency FundingCcy   = *static_cast< ARM_ZeroCurve* >(GetMktDataManager()->GetData(GetKeys()[FundingKey]))->GetCurrencyUnit();
    SetFundingCcy(FundingCcy);

    ARM_Currency BasisCcy     = *static_cast< ARM_ZeroCurve* >(GetMktDataManager()->GetData(GetKeys()[BasisKey]))->GetCurrencyUnit();
    SetBasisCcy(BasisCcy);

    if(string(cpnCcy->GetCcyName()) == string(FundingCcy.GetCcyName()))
    {
        // no basis, no forex
        SetForeignCcy(*cpnCcy);
        SetDomesticCcy(*cpnCcy);
    }
    else
    {
        ARM_Forex* forex = static_cast< ARM_Forex* >(GetMktDataManager()->GetData(GetKeys()[ForexKey]));
        SetForeignCcy(*forex->GetMainCurrency());
        SetDomesticCcy(*forex->GetMoneyCurrency());
    }

	/// Check input datas
	CheckDataAndTimeIt();

	if (itsCalibND == TWODCalib)
		itsMrsCalib = false;

	/// Create the Generic Security paid in coupon currency
	CreateAndSetDealDescriptionAndTimeIt();

	ARM_GP_Vector* resetDates = itsCalibDateStrip->GetResetDates();
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	int resetSize=resetDates->size();
	int i=0;
	while( (i < resetSize) && (const_cast< ARM_Curve& >(GetFixCpn()).Interpolate((*resetDates)[i]-asOfDate) > 0.0) )
		++i;
	SetNbFixFlows(i);

	/// Create the SFRM pricing model with its default parameters
	CreateAndSetModelAndTimeIt();

	/// Create the calibration set for volatility bootstapping and
	/// strike spread adjustments
	CreateAndSetCalibrationAndTimeIt();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CSOCalculator
///	Routine: Destructor
///	Returns: void
///	Action : destroys the object
/////////////////////////////////////////////////////////////////
ARM_CSOCalculator::~ARM_CSOCalculator()
{	
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_CSOCalculator
///	Routine: Copy constructor
///	Returns: a built object
///	Action : builds the object
/////////////////////////////////////////////////////////////////
ARM_CSOCalculator::ARM_CSOCalculator(const ARM_CSOCalculator& rhs)
: ARM_GenCSOCalculator(rhs),
	itsSigmaCalib(rhs.itsSigmaCalib),
	itsCalibND(rhs.itsCalibND),
	itsMrsCalib(rhs.itsMrsCalib),
	itsThetaCalib(rhs.itsThetaCalib),
	itsModelDatas(rhs.itsModelDatas),
	itsUpdateInitCalib(rhs.itsUpdateInitCalib)
{
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CSOCalculator
///	Routine: assignment operator
///	Returns: 
///	Action : 
/////////////////////////////////////////////////////////////////
ARM_CSOCalculator& ARM_CSOCalculator::operator=(const ARM_CSOCalculator& rhs)
{
	if (this != &rhs)
	{
		ARM_GenCalculator::operator=(rhs);

		itsSigmaCalib				= rhs.itsSigmaCalib;
		itsCalibND					= rhs.itsCalibND;
		itsMrsCalib					= rhs.itsMrsCalib;
		itsThetaCalib				= rhs.itsThetaCalib;
		itsModelDatas				= rhs.itsModelDatas;
		itsUpdateInitCalib			= rhs.itsUpdateInitCalib;
	}

	return *this;
}

////////////////////////////////////////////////////
///	Class   : ARM_CSOCalculator
///	Routines: Clone,View
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_CSOCalculator::Clone() const
{
	return new ARM_CSOCalculator(*this);
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: InitPriceableColumns
///	Returns: nothing
///	Action : Feed the MDM with input market objects
/////////////////////////////////////////////////////////////////
void ARM_CSOCalculator::InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const
{
    string zeroValue("0");

	rowDescVec[CSO] = zeroValue;
    rowTypeVec[CSO] = ARM_DOUBLE;
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_CSOCalculator
///	Routine: MiddleRows
///	Returns: ARM_RowInfo
///	Action : create a row of a deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_CSOCalculator::MiddleRows( size_t eventIdx, const ARM_DateStripCombiner& datesStructure ) const
{
    ARM_DateStripPtr exerDateStrip = datesStructure.GetDateStrip(EXERCISE_SCHED);
	ARM_DateStripPtr cpnDateStrip = itsStructDateStrip;
	ARM_DateStripPtr fundDateStrip = itsStructDateStrip;

    size_t eventSize = exerDateStrip->GetResetDates()->size();
    size_t descSize = sizeof(CSOColNamesTable)/sizeof(CSOColNamesTable[0]);

	bool isFirstEvent = (eventIdx==itsFirstEventIdx);

    string modelName = GetKeys()[YcKey];

	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();   

    vector< string > rowDescVec(descSize);
    vector< ARM_GP_VALUE_TYPE > rowTypeVec(descSize, ARM_MISSING_TYPE); 

    /// Set default 0 value for each column to be able to sum it
    InitPriceableColumns(rowDescVec,rowTypeVec);

	// EventDate 
    double eventDate=(*(exerDateStrip->GetResetDates()))[eventIdx];
    CC_Ostringstream eventDateDesc;
    eventDateDesc << CC_NS(std,fixed) << eventDate;
    rowDescVec[EventDate] = eventDateDesc.str();
    rowTypeVec[EventDate] = ARM_DATE_TYPE;

    double startDate=(*(itsStructDateStrip->GetFlowStartDates()))[eventIdx];
    CC_Ostringstream startDateDesc;
    startDateDesc << CC_NS(std,fixed) << startDate;
    rowDescVec[StartDate] = startDateDesc.str();
    rowTypeVec[StartDate] = ARM_DATE_TYPE;

    double payDate=(*(itsStructDateStrip->GetPaymentDates()))[eventIdx];
    CC_Ostringstream payDateDesc;
    payDateDesc << CC_NS(std,fixed) << payDate;
    rowDescVec[PayDate] = payDateDesc.str();
    rowTypeVec[PayDate] = ARM_DATE_TYPE;

	ARM_Date cmsDate1(startDate), cmsDate2(startDate);

	if (IsCMSIndex((ARM_INDEX_TYPE)GetCMSLong()))
		cmsDate1.AddMonths( (GetCMSLong()-K_CMS1+1)*12 );
	else
	{
		int tmpFreq = FromLiborTypeToFrequency(GetCMSLong());
		cmsDate1.AddMonths( 12 / tmpFreq );
	}

	if (IsCMSIndex((ARM_INDEX_TYPE)GetCMSShort()))
		cmsDate2.AddMonths( (GetCMSShort()-K_CMS1+1)*12 );
	else
	{
		int tmpFreq = FromLiborTypeToFrequency(GetCMSShort());
		cmsDate2.AddMonths( 12 / tmpFreq );
	}

    CC_Ostringstream endDateCMS1Desc;
    endDateCMS1Desc << CC_NS(std,fixed) << cmsDate1.GetJulian();
    rowDescVec[EndDateCMS1] = endDateCMS1Desc.str();
    rowTypeVec[EndDateCMS1] = ARM_DATE_TYPE;

    CC_Ostringstream endDateCMS2Desc;
    endDateCMS2Desc << CC_NS(std,fixed) << cmsDate2.GetJulian();
    rowDescVec[EndDateCMS2] = endDateCMS2Desc.str();
    rowTypeVec[EndDateCMS2] = ARM_DATE_TYPE;

	//IT
	CC_Ostringstream ITDesc;
	ITDesc << CC_NS(std,fixed) << (*(itsStructDateStrip->GetInterestTerms()))[eventIdx];
	rowDescVec[IT] = ITDesc.str();
	rowTypeVec[IT] = ARM_DOUBLE;

	//CMS1
	CC_Ostringstream CMS1Desc;
	CMS1Desc << "SWAPRATE(" << modelName << "," << CSOColNamesTable[StartDate] << "[i],";
	CMS1Desc << CSOColNamesTable[EndDateCMS1] << "[i])";
	rowDescVec[CMS1] = CMS1Desc.str();
	rowTypeVec[CMS1] = ARM_STRING;

	//CMS2
	CC_Ostringstream CMS2Desc;
	CMS2Desc << "SWAPRATE(" << modelName << "," << CSOColNamesTable[StartDate] << "[i],";
	CMS2Desc << CSOColNamesTable[EndDateCMS2] << "[i])";
	rowDescVec[CMS2] = CMS2Desc.str();
	rowTypeVec[CMS2] = ARM_STRING;

	double fixCpn = const_cast< ARM_Curve& >(GetFixCpn()).Interpolate(eventDate-asOfDate);

	//CpnFix
	CC_Ostringstream cpnFixDesc;
	cpnFixDesc << CC_NS(std,fixed) << fixCpn;
	rowDescVec[CpnFix] = cpnFixDesc.str();
	rowTypeVec[CpnFix] = ARM_DOUBLE;

	//CpnMin
	CC_Ostringstream strikeDesc;
	strikeDesc << CC_NS(std,fixed) << const_cast< ARM_Curve& >(GetCpnFloor()).Interpolate(eventDate-asOfDate);
	rowDescVec[Strike] = strikeDesc.str();
	rowTypeVec[Strike] = ARM_DOUBLE;

	//Notional
	CC_Ostringstream notDesc;
	notDesc << CC_NS(std,fixed) << const_cast< ARM_Curve& >(GetCpnNominal()).Interpolate(eventDate-asOfDate);
	rowDescVec[Notional] = notDesc.str();
	rowTypeVec[Notional] = ARM_DOUBLE;

	//Leverage
	CC_Ostringstream leverageDesc;
	leverageDesc << CC_NS(std,fixed) << const_cast< ARM_Curve& >(GetLeverageLong()).Interpolate(eventDate-asOfDate);
	rowDescVec[Leverage] = leverageDesc.str();
	rowTypeVec[Leverage] = ARM_DOUBLE;

	//spreadoptionlet
	CC_Ostringstream spreadoptionletDesc;
	if (fixCpn != 0.0)
	{
		spreadoptionletDesc << CSOColNamesTable[CpnFix] << "[i]*" << CSOColNamesTable[IT] << "[i]";
		spreadoptionletDesc << "*DF(" << modelName << "," << CSOColNamesTable[PayDate] << "[i])*";
		spreadoptionletDesc << CSOColNamesTable[Notional] << "[i]";
	}
	else
	{
		spreadoptionletDesc << "(MAX(" << CSOColNamesTable[Leverage] << "[i]";
		spreadoptionletDesc << "*(" << CSOColNamesTable[CMS1] << "[i]-" << CSOColNamesTable[CMS2] << "[i])-" << CSOColNamesTable[Strike] << "[i],0)+";
		spreadoptionletDesc << CSOColNamesTable[Strike] << "[i])*" << CSOColNamesTable[IT] << "[i]*DF(" << modelName << ",";
		spreadoptionletDesc << CSOColNamesTable[PayDate] << "[i])*" << CSOColNamesTable[Notional] << "[i]";
	}
	rowDescVec[SpreadOptionLet] = spreadoptionletDesc.str();
	rowTypeVec[SpreadOptionLet] = ARM_STRING;

	bool isLastEvent = (eventIdx==(eventSize-1));
	CC_Ostringstream spreadoptionDesc;

	if (!isLastEvent)
	{
		spreadoptionDesc << CSOColNamesTable[SpreadOptionLet] << "[i]+PV(";
		spreadoptionDesc << CSOColNamesTable[SpreadOption] << "[i+1])";
	}
	else
	{
		spreadoptionDesc << CSOColNamesTable[SpreadOptionLet] << "[i]";
	}
	rowDescVec[SpreadOption] = spreadoptionDesc.str();
	rowTypeVec[SpreadOption] = ARM_STRING;

	CC_Ostringstream strikeflowDesc;
	strikeflowDesc << CSOColNamesTable[Strike] << "[i]*" << CSOColNamesTable[IT] << "[i]*DF(" << modelName;
	strikeflowDesc << "," << CSOColNamesTable[PayDate] << "[i])*" << CSOColNamesTable[Notional] << "[i]";
	rowDescVec[StrikeFlow] = strikeflowDesc.str();
	rowTypeVec[StrikeFlow] = ARM_STRING;

	CC_Ostringstream pvstrikeflowDesc;

	if (!isLastEvent)
	{
		pvstrikeflowDesc << CSOColNamesTable[StrikeFlow] << "[i]+PV(";
		pvstrikeflowDesc << CSOColNamesTable[PVStrikeFlow] << "[i+1])";
	}
	else
	{
		pvstrikeflowDesc << CSOColNamesTable[StrikeFlow] << "[i]";
	}
	rowDescVec[PVStrikeFlow] = pvstrikeflowDesc.str();
	rowTypeVec[PVStrikeFlow] = ARM_STRING;

    double fundStartDate=(*(fundDateStrip->GetFlowStartDates()))[eventIdx];
    CC_Ostringstream fundStartDateDesc;
    fundStartDateDesc << CC_NS(std,fixed) << fundStartDate;
    rowDescVec[FundStartDate] = fundStartDateDesc.str();
    rowTypeVec[FundStartDate] = ARM_DATE_TYPE;

    double fundEndDate=(*(fundDateStrip->GetFlowEndDates()))[eventIdx];
    CC_Ostringstream fundEndDateDesc;
    fundEndDateDesc << CC_NS(std,fixed) << fundEndDate;
    rowDescVec[FundEndDate] = fundEndDateDesc.str();
    rowTypeVec[FundEndDate] = ARM_DATE_TYPE;

	//fundSpread
	CC_Ostringstream fundspreadDesc;
	fundspreadDesc << CC_NS(std,fixed) << const_cast< ARM_Curve& >(GetFundSpread()).Interpolate(eventDate-asOfDate);
	rowDescVec[FundSpread] = fundspreadDesc.str();
	rowTypeVec[FundSpread] = ARM_DOUBLE;

	string stringMatu(YearTermToStringMatu(1./GetFundFreq()));

	CC_Ostringstream fundFragmentNPVDesc;
	fundFragmentNPVDesc << "SWAP(" << modelName << "," << CSOColNamesTable[FundStartDate] << "[i]," << CSOColNamesTable[FundEndDate] << "[i],0.0,P,,,";
	fundFragmentNPVDesc << stringMatu << ",," << CSOColNamesTable[FundSpread] << "[i])*" << CSOColNamesTable[Notional] << "[i]";
	rowDescVec[FundFragmentNPV] = fundFragmentNPVDesc.str();
	rowTypeVec[FundFragmentNPV] = ARM_STRING;

	CC_Ostringstream fundnpvDesc;

	if (!isLastEvent)
	{
		fundnpvDesc << CSOColNamesTable[FundFragmentNPV] << "[i]+PV(";
		fundnpvDesc << CSOColNamesTable[FundNPV] << "[i+1])";
	}
	else
	{
		fundnpvDesc << CSOColNamesTable[FundFragmentNPV] << "[i]";
	}
	rowDescVec[FundNPV] = fundnpvDesc.str();
	rowTypeVec[FundNPV] = ARM_STRING;

	//spreadoptionlet
	CC_Ostringstream flowDesc;
	if (fixCpn != 0.0)
	{
		flowDesc << CSOColNamesTable[CpnFix] << "[i]*DF(" << modelName << "," << CSOColNamesTable[PayDate] << "[i])*";
		flowDesc << CSOColNamesTable[Notional] << "[i]";
	}
	else
	{
		flowDesc << CSOColNamesTable[SpreadOptionLet] << "[i]";
	}
	rowDescVec[Flow] = flowDesc.str();
	rowTypeVec[Flow] = ARM_STRING;

	CC_Ostringstream totalDesc;

	if (!isLastEvent)
	{
		totalDesc << CSOColNamesTable[Flow] << "[i]+PV(";
		totalDesc << CSOColNamesTable[TotalFlow] << "[i+1])";
	}
	else
	{
		totalDesc << CSOColNamesTable[Flow] << "[i]";
	}
	rowDescVec[TotalFlow] = totalDesc.str();
	rowTypeVec[TotalFlow] = ARM_STRING;

	//fees
	CC_Ostringstream feesDesc;
	feesDesc << CC_NS(std,fixed) << const_cast< ARM_Curve& >(GetFees()).Interpolate(eventDate-asOfDate);
	rowDescVec[Fees] = feesDesc.str();
	rowTypeVec[Fees] = ARM_DOUBLE;

	CC_Ostringstream bermudaProfileDesc;

	if (!isLastEvent)
	{
		bermudaProfileDesc << "MAX(" << CSOColNamesTable[TotalFlow] << "[i]-" << CSOColNamesTable[FundNPV] << "[i]-" << CSOColNamesTable[Fees] << "[i],";
		bermudaProfileDesc << "PV(" << CSOColNamesTable[BermudaProfile] << "[i+1]))";
	}
	else
	{
		bermudaProfileDesc << "MAX(" << CSOColNamesTable[TotalFlow] << "[i]-" << CSOColNamesTable[FundNPV] << "[i]-" << CSOColNamesTable[Fees] << "[i],0)";
	}
	rowDescVec[BermudaProfile] = bermudaProfileDesc.str();
	rowTypeVec[BermudaProfile] = ARM_STRING;

	CC_Ostringstream csoDesc;
	if (isFirstEvent)
	{
		csoDesc << "MAX(" << CSOColNamesTable[TotalFlow] << "[i]-" << CSOColNamesTable[FundNPV] << "[i]-" << CSOColNamesTable[Fees] << "[i],";
		csoDesc << "PV(" << CSOColNamesTable[BermudaProfile] << "[i+1]))";
		rowDescVec[CSO] = csoDesc.str();
		rowTypeVec[CSO] = ARM_STRING;
	}

	return ARM_RowInfo(rowDescVec,rowTypeVec);
}




/////////////////////////////////////////////////////////////////
///	Class  : ARM_CSOCalculator
///	Routine: CreateAndSetModel
///	Returns: void
///	Action : creates the model
/////////////////////////////////////////////////////////////////
void ARM_CSOCalculator::CreateAndSetModel()
{
    /// Get yield curve
	ARM_ZeroCurve* curve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));

	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
//	const ARM_DateStripCombiner& dateStructure = DatesStructure();

	ARM_INDEX_TYPE liborType = GetIndexType();
	ARM_IRIndex IRIndex(liborType,K_DEF_FREQ,K_DEF_FREQ,GetCurrencyUnit());

	/// Creates default values for volatility & mean reversion
    /// (calibration will be called at pricing time) and set them
    /// in the reference model
	ARM_GP_Vector* resetDates = itsCalibDateStrip->GetResetDates();
	double firstCalibTime = (*resetDates)[GetNbFixFlows()]-asOfDate;
	 
	ARM_GP_Vector defaultTimes(1,firstCalibTime);
	ARM_GP_Vector defaultSigmas(1,SIGMA_DEFAULT_VALUE);
	ARM_GP_Vector defaultBetas(1,BETA_DEFAULT_VALUE);

	ARM_CurveModelParam volParam( ARM_ModelParamType::Volatility, &defaultSigmas, &defaultTimes );
	ARM_CurveModelParam newBetaParam( ARM_ModelParamType::Beta, &defaultBetas, &defaultTimes );
	ARM_CurveModelParam* mrsParam = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(GetKeys()[MrsKey]));
    if(!mrsParam || mrsParam->GetType() != ARM_ModelParamType::MeanReversion)
    {
		CC_Ostringstream os;
		os << ARM_USERNAME << " : MRS Param for key=" << GetKeys()[MrsKey] << " is expected in the Market Data Manager";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
    }

	ARM_ModelParam* correlParam = dynamic_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[CorrelKey]));
    if(!correlParam || correlParam->GetType() != ARM_ModelParamType::Correlation)
    {
		CC_Ostringstream os;
		os << ARM_USERNAME << " : Correl Param for key=" << GetKeys()[CorrelKey] << " is expected in the Market Data Manager";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
    }

	ARM_ModelParamVector paramVector(4);
	paramVector[0] = &volParam;
	paramVector[1] = mrsParam;
	paramVector[2] = correlParam;

	/// use the factory class to get the model params!
	ARM_ModelParamsSFRM* pSFRMModelParams = 
		ARM_ModelParamsSFRMFactory.Instance()->CreateModelParamsSFRM(paramVector,&IRIndex,SFRM_NB_FACTORS,SFRM_VOL_TYPE);

	/// Build the default stochastic model of the calculator : SFRM 2F
    ARM_PricingModelPtr refModel( ARM_PricingModelPtr(static_cast< ARM_PricingModel* >(new ARM_SFRM( CreateClonedPtr( curve ), *pSFRMModelParams))) );

	// Delte the model params because it is cloned in the model
	delete pSFRMModelParams;

    /// Create a Numeraire and set it
    ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::TerminalZc ) );
    refModel->SetNumeraire(numeraire);

	const ARM_DealDescription dealDescr = GetGenSecurity()->GetDealDescription();
    double lastEventTime = atof(dealDescr.GetElem(dealDescr.GetRowsNb()-1,EventDate).c_str()) - asOfDate;
    int nbSteps=static_cast<int>(floor(TREE_NBSTEPS_PER_YEAR*lastEventTime/K_YEAR_LEN));
	int initialStep =static_cast<int> (nbSteps/10.0);

    int schedulerType=ARM_SchedulerBase::ConstantVariance;
    ARM_GP_Vector schedulerDatas(3);
    schedulerDatas[0] = nbSteps;
    schedulerDatas[1] = initialStep;
    int samplerType=ARM_SamplerBase::NormalCentred;
    ARM_GP_Vector samplerDatas;
    int truncatorType=ARM_TruncatorBase::StandardDeviation;
    ARM_GP_Vector truncatorDatas(1,STD_DEV_RATIO);
    int reconnectorType=ARM_ReconnectorBase::Mean;
    int smootherType=ARM_SmootherBase::DoNothing;
    ARM_TreeBase* tree = ARM_TreeFactory.Instance()->CreateTreeND(SFRM_NB_FACTORS,schedulerType,schedulerDatas,
        samplerType,samplerDatas,truncatorType,truncatorDatas,false,reconnectorType,smootherType);
    refModel->SetNumMethod(ARM_NumMethodPtr( tree ) );

	/// Set the model
	SetPricingModel(refModel);

}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CSOCalculator
///	Routine: ColumnNames
///	Returns: ARM_RowInfo
///	Action : create the column names of the deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_CSOCalculator::ColumnNames() const
{
	size_t colNamesSize = sizeof(CSOColNamesTable)/sizeof(CSOColNamesTable[0]);
    vector< string > colNamesVec(colNamesSize);
    vector< ARM_GP_VALUE_TYPE > colTypeVec(colNamesSize, ARM_STRING); 

    for(size_t i=0;i<colNamesSize; ++i)
        colNamesVec[i] = CSOColNamesTable[i];

    ARM_RowInfo rowInfo(colNamesVec,colTypeVec);

    return rowInfo;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CSOCalculator
///	Routine: ComputePricingData
///	Returns: a ARM_MultiTypeDict
///	Action : get pricing data of the calculaltor
/////////////////////////////////////////////////////////////////

void ARM_CSOCalculator::ComputePricingData() const
{
	if (!itsHasBeenPriced)
		const_cast<ARM_CSOCalculator*>(this)->PriceAndTimeIt();

	if (itsProductsToPrice[CSOPrice])
		GetPricingData()[ "CSOPrice" ] = itsCSOPrice;
	
	if (itsProductsToPrice[structPrice])
		GetPricingData()[ "SOPrice" ] = itsStructPrice;		

}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CSOCalculator
///	Routine: GenerateEquivalentSwoptStrikes
///	Returns: ARM_GP_Vector*
///	Action : generate equivalent swaptionStrikes on fund 
/////////////////////////////////////////////////////////////////
void ARM_CSOCalculator::ComputeEquivalentSwoptStrikes(ARM_GP_Vector& equivStrikes1, ARM_GP_Vector& equivStrikes2)
{
	ARM_GP_Vector* startDates = itsCalibDateStrip->GetFlowStartDates();
	ARM_GP_Vector* resetDates = itsCalibDateStrip->GetResetDates();
	int fixFreq = GetCurrencyUnit()->GetFixedPayFreq();
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

	ARM_ZeroCurve* zc = dynamic_cast< ARM_ZeroCurve* >( GetMktDataManager()->GetData(GetKeys()[YcKey]) );
	ARM_BSModel* SOBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[SoModelKey]) );

	list< ARM_Security* > spreadOptionList;

	int size = startDates->size() - GetNbFixFlows();

	for (int i = 0; i < size; i++)
	{
		double fixCpn = const_cast< ARM_Curve& >(GetCpnFloor()).Interpolate((*resetDates)[i+GetNbFixFlows()]-asOfDate);
		double leverage = const_cast< ARM_Curve& >(GetLeverageLong()).Interpolate((*resetDates)[i+GetNbFixFlows()]-asOfDate);
		double strike = fixCpn/leverage*100.;

		ARM_Date startDate((*startDates)[i+GetNbFixFlows()]);
		ARM_Date endDate((*startDates)[i+GetNbFixFlows()]);
		endDate.AddMonths(12/GetCpnFreq());

		ARM_Vector* fixing1 = NULL;
		ARM_Vector* fixing2 = NULL;

		ARM_ReferenceValue	vWeight1(1.0);
		ARM_ReferenceValue	vWeight2(1.0);
		ARM_ReferenceValue  aStrike(strike);

		ARM_SpreadOption so(startDate,
						   endDate,
						   K_CAP,
						   &aStrike,
						   (ARM_INDEX_TYPE)GetCMSShort(),
						   (ARM_INDEX_TYPE)GetCMSLong(),
						   &vWeight1, &vWeight2,
						   GetCpnDaycount(),
						   GetCpnFreq(),
						   GetCpnFreq(),
						   K_ADVANCE,
						   K_ARREARS,
						   GetCurrencyUnit(),
						   10000,
						   fixing1,
						   fixing2);

		so.SetModel(SOBSModel);

		spreadOptionList.push_back(static_cast<ARM_SpreadOption*>(so.Clone()));

		double firstFwd = so.GetSpreadLeg()->GetFirstLeg()->GetFwdRates()->Elt(0);
		double secondFwd = so.GetSpreadLeg()->GetSecondLeg()->GetFwdRates()->Elt(0);

		double firstVol = so.GetSpreadLeg()->GetFirstLeg()->GetVolatilityFwdRates()->Elt(0);
		double secondVol = so.GetSpreadLeg()->GetSecondLeg()->GetVolatilityFwdRates()->Elt(0);

		equivStrikes2[i] = firstFwd + so.CptStrike(firstFwd,secondFwd,1.0,1.0,firstVol,secondVol,strike,1);
		equivStrikes1[i] = secondFwd + so.CptStrike(firstFwd,secondFwd,1.0,1.0,firstVol,secondVol,strike,2);
	}

	itsSpreadOptionPF = ARM_StdPortfolioPtr(new ARM_StdPortfolio(spreadOptionList));
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CSOCalculator
///	Routine: CreateCSOSwaptionCap
///	Returns: void
///	Action : create the swaption PF
/////////////////////////////////////////////////////////////////
void ARM_CSOCalculator::CreateCSOSwaptionCap()
{
    ARM_GP_Vector* resetDates = itsCalibDateStrip->GetResetDates();
	ARM_GP_Vector* startDates = itsCalibDateStrip->GetFlowStartDates();

	int resetSize = resetDates->size()-GetNbFixFlows();

	ARM_GP_Vector equivStrikes1(resetSize);
	ARM_GP_Vector equivStrikes2(resetSize);
	ComputeEquivalentSwoptStrikes(equivStrikes1,equivStrikes2);

	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

	int portfolioSize = resetSize;

	//ARM_GP_Vector initBreakpointTimes(portfolioSize);
	//ARM_GP_Vector initVolatility(portfolioSize,SIGMA_DEFAULT_VALUE);	

	list< ARM_Security* > swaption1List;
	list< ARM_Security* > swaption2List;
	

	for (int i = 0; i < portfolioSize; ++i)
	{
		double resetDate = (*resetDates)[i+GetNbFixFlows()];
		ARM_Date startDate((*startDates)[i+GetNbFixFlows()]);
		ARM_Date expiryDate(resetDate);
		//initBreakpointTimes[i]=resetDate-asOfDate;

		ARM_Date endDate1(startDate);

		if (IsCMSIndex((ARM_INDEX_TYPE)GetCMSLong()))
			endDate1.AddMonths( (GetCMSLong()-K_CMS1+1)*12 );
		else
		{
			int tmpFreq = FromLiborTypeToFrequency(GetCMSLong());
			endDate1.AddMonths( 12 / tmpFreq );
		}

		ARM_Swap stdSwap1(startDate,
						 endDate1,
						 GetIndexType(),0.0,1.0,K_RCV,K_DEF_FREQ,K_DEF_FREQ,GetCurrencyUnit());

		double equivStrike1 = equivStrikes1[i];

		int RecOrPay = K_RCV;
		ARM_Swaption swaption1((&stdSwap1),RecOrPay,K_EUROPEAN,equivStrike1,expiryDate);
		
		swaption1List.push_back(static_cast<ARM_Swaption*>(swaption1.Clone()));

		ARM_Date endDate2(startDate);

		if (IsCMSIndex((ARM_INDEX_TYPE)GetCMSShort()))
			endDate2.AddMonths( (GetCMSShort()-K_CMS1+1)*12 );
		else
		{
			int tmpFreq = FromLiborTypeToFrequency(GetCMSShort());
			endDate2.AddMonths( 12 / tmpFreq );
		}

		ARM_Swap stdSwap2(startDate,
						 endDate2,
						 GetIndexType(),0.0,1.0,K_RCV,K_DEF_FREQ,K_DEF_FREQ,GetCurrencyUnit());

		double equivStrike2 = equivStrikes2[i];

		ARM_Swaption swaption2((&stdSwap2),RecOrPay,K_EUROPEAN,equivStrike2,expiryDate);
		
		swaption2List.push_back(static_cast<ARM_Swaption*>(swaption2.Clone()));
	}

	itsSwaptionCMS1PF = ARM_StdPortfolioPtr(new ARM_StdPortfolio(swaption1List));
	itsSwaptionCMS2PF = ARM_StdPortfolioPtr(new ARM_StdPortfolio(swaption2List));
	
	for(i=0;i<portfolioSize;++i)
	{
        itsSwaptionCMS1PF->SetWeight(DEFAULT_WEIGHT,i);
		itsSwaptionCMS1PF->SetPrice(DEFAULT_PRICE,i);
		itsSwaptionCMS1PF->SetPrecision(DEFAULT_PRECISION,i);
		
        itsSwaptionCMS2PF->SetWeight(DEFAULT_WEIGHT,i);
		itsSwaptionCMS2PF->SetPrice(DEFAULT_PRICE,i);
		itsSwaptionCMS2PF->SetPrecision(DEFAULT_PRECISION,i);
		
		list< ARM_Security* > swaption1_2List;
		swaption1_2List.push_back(static_cast<ARM_Swaption*>((*itsSwaptionCMS1PF).GetAsset(i)->Clone()));
		swaption1_2List.push_back(static_cast<ARM_Swaption*>((*itsSwaptionCMS2PF).GetAsset(i)->Clone()));

		(itsSwaptionCMSPF_1_2).push_back(ARM_StdPortfolioPtr(new ARM_StdPortfolio(swaption1_2List)));
	}

	/*ARM_SFRM* refModel;
	refModel = dynamic_cast< ARM_SFRM* >(&*GetPricingModel());

	((ARM_CurveModelParam&) refModel->GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility)).SetValuesAndTimes(&initBreakpointTimes,&initVolatility); 
*/
}




/////////////////////////////////////////////////////////////////
///	Class  : ARM_CSOCalculator
///	Routine: GetCFCalibMethod
///	Returns: ARM_CalibMethod
///	Action : get calibration method for sigma
/////////////////////////////////////////////////////////////////
ARM_CalibMethod* ARM_CSOCalculator::GetSigmaCalibMethod() const   
{
	if ( (itsCalibND == ONEDCalib) || (itsMrsCalib == true) )
	{
		if(itsMrsCalib == true)
		{
			return GetCalibMethod()->GetlinkedMethod()->GetlinkedMethod();
		}
		else
		{
			return GetCalibMethod()->GetlinkedMethod();
		}
	}
	else
	{
		return GetCalibMethod()->GetlinkedMethod();
	}
}




/////////////////////////////////////////////////////////////////
///	Class  : ARM_CSOCalculator
///	Routine: CreateDiagonalSwaption
///	Returns: a portfolio
///	Action : create the list of diagonal swaptions
/////////////////////////////////////////////////////////////////
void ARM_CSOCalculator::CreateDiagonalSwaption()
{
	ARM_GP_Vector* resetDates = itsCalibDateStrip->GetResetDates();
	ARM_GP_Vector* startDates = itsCalibDateStrip->GetFlowStartDates();

	ARM_BSModel* CFBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[CfModelKey]) );

	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

	int portfolioSize = resetDates->size() - GetNbFixFlows();

	int i;

	list< ARM_Security* > swaptionList;
	
    ARM_Swaption swaption;

	for (i = 0; i < portfolioSize; ++i)
	{
		double resetDate = (*resetDates)[i+GetNbFixFlows()];
		ARM_Date startDate((*startDates)[i+GetNbFixFlows()]);
		ARM_Date expiryDate(resetDate);

		ARM_Swap stdSwap(startDate,
                GetEndDate(),
                GetIndexType(),0.0,1.0,K_RCV,K_DEF_FREQ,K_DEF_FREQ,GetCurrencyUnit());

		stdSwap.SetModel(CFBSModel);
        
		double equivStrike = stdSwap.PriceToRate((ARM_Date) asOfDate, 0.0);

		int RecOrPay = K_RCV;
		ARM_Swaption swaption((&stdSwap),RecOrPay,K_EUROPEAN,equivStrike,expiryDate);
		
		swaptionList.push_back(static_cast<ARM_Swaption*>(swaption.Clone()));
	}

	itsSwaptionPF = ARM_StdPortfolioPtr(new ARM_StdPortfolio(swaptionList));
    
	
	for(i=0;i<itsSwaptionPF->size();++i)
	{
		itsSwaptionPF->SetPrecision(MRS_DEFAULT_PRECISION,i);
        itsSwaptionPF->SetWeight(DEFAULT_WEIGHT,i);
		itsSwaptionPF->SetPrice(DEFAULT_PRICE,i);
	}
}

////////////////////////////////////////////////////////////////
///	Class  : ARM_CSOCalculator
///	Routine: computeDIAGSwaptionPortfolioPrices
///	Returns: nothing
///	Action : compute market target prices of the swaption portfolio
/////////////////////////////////////////////////////////////////
void ARM_CSOCalculator::computeSOPortfolioPrices()
{
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	
	ARM_BSModel* SOBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[SoModelKey]) );

	ARM_Date cmsDate1(GetEndDate());
	if (IsCMSIndex((ARM_INDEX_TYPE)GetCMSLong()))
		cmsDate1.AddMonths( (GetCMSLong()-K_CMS1+1)*12 );

	ARM_DateStrip datestrip(GetStartDate(),cmsDate1,CORREL_FREQ,GetCpnDaycount(),
						GETDEFAULTVALUESTR,
						K_MOD_FOLLOWING,K_ADJUSTED,K_SHORTSTART,
						GETDEFAULTVALUE,CORREL_FREQ);

	ARM_GP_Vector lowerBound(1,CORREL_LOWER_BOUND);
	ARM_GP_Vector upperBound(1,CORREL_UPPER_BOUND);

	ARM_CorrelTrigoMatParam* correlparam = new ARM_CorrelTrigoMatParam(CORREL_DEFAULT_INITVALUE,
															datestrip,
															asOfDate,
															"LINEAR",
															&lowerBound,
															&upperBound);

	ARM_StdPortfolioPtr portfolioSO = GetCalibMethod()->GetPortfolio();
	ARM_SpreadOption* so;

	int sizeSO = (*portfolioSO).size();

	for (int i = 0; i < sizeSO; i++)
	{
		so = static_cast< ARM_SpreadOption* >(portfolioSO->GetAsset(i));
		so->SetModel(SOBSModel);
		double SO_price = so->ComputePrice();
		portfolioSO->SetPrice(SO_price,i);
		portfolioSO->SetWeight(DEFAULT_WEIGHT,i);
		portfolioSO->SetPrecision(DEFAULT_PRECISION, i);
	}

	int paramSize=GetCalibMethod()->GetCalibParams().size();

	if(paramSize == 0 )
		GetCalibMethod()->GetCalibParams().push_back(correlparam);
	else
	{
		delete GetCalibMethod()->GetCalibParam(0);
		(GetCalibMethod()->GetCalibParams())[0]= correlparam;
	}
}
	

////////////////////////////////////////////////////////////////
///	Class  : ARM_CSOCalculator
///	Routine: computeDIAGSwaptionPortfolioPrices
///	Returns: nothing
///	Action : compute market target prices of the swaption portfolio
/////////////////////////////////////////////////////////////////
void ARM_CSOCalculator::computeDIAGSwaptionPortfolioPrices()
{
    double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	ARM_BSModel* oswBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswModelKey]) );

	ARM_CalibMethod* mrsCalibMethod =(GetCalibMethod())->GetlinkedMethod();
    ARM_StdPortfolioPtr swaptionPF = mrsCalibMethod->GetPortfolio();
    int nbOSW=itsSwaptionPF->GetSize();

	ARM_Swaption* swaption;
	double price;
	
	for(int i=0;i<nbOSW;++i)
	{
		swaption = static_cast< ARM_Swaption* >(swaptionPF->GetAsset(i));
		swaption->SetModel(oswBSModel);
		price=swaption->ComputePrice();
		swaptionPF->SetPrice(price,i);
	}

	size_t mrsIdx,paramSize=mrsCalibMethod->GetCalibParams().size();
    for(mrsIdx=0;mrsIdx<paramSize;++mrsIdx)
    {
        if( mrsCalibMethod->GetCalibParam(mrsIdx) &&
            (mrsCalibMethod->GetCalibParams())[mrsIdx]->GetType() == ARM_ModelParamType::MeanReversion )
            break;
    }        

	ARM_GP_Vector initTimes(1,0.0);
	ARM_GP_Vector initMRS(1,MRS_DEFAULT_VALUE);
	ARM_GP_Vector mrsLowerBound(1,MRS_LOWER_BOUND);
	ARM_GP_Vector mrsUpperBound(1,MRS_UPPER_BOUND);

	ARM_CurveModelParam* mrs = new ARM_CurveModelParam(
		ARM_ModelParamType::MeanReversion,
		&initMRS,
		&initTimes,
		"MEANREVERSION",
		"STEPUPRIGHT",
		&mrsLowerBound,
		&mrsUpperBound);

	if(mrsIdx >= paramSize || paramSize == 0)
		mrsCalibMethod->GetCalibParams().push_back(mrs);
    else
    {
        delete mrsCalibMethod->GetCalibParam(mrsIdx);
        (mrsCalibMethod->GetCalibParams())[mrsIdx] = mrs;
    }
}




/////////////////////////////////////////////////////////////////
///	Class  : ARM_CSOCalculator
///	Routine: CreateEmptyCalibration
///	Returns: void
///	Action : create an empty calibration structure
/////////////////////////////////////////////////////////////////
ARM_CalibMethodPtr ARM_CSOCalculator::CreateSigmaCalibration()
	
{
		CreateCSOSwaptionCap();
		ARM_CalibMethod* volCalibMethod = NULL;
		if (itsSigmaCalib == BootstrapCalib)
		{
			volCalibMethod = new ARM_CalibMethod(itsSwaptionCMS1PF,ARM_ModelParamVector(),ARM_CalibMethodType::Bootstrap1D);
		}
		else
		{
			double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
			ARM_GP_Vector* resetDates = itsCalibDateStrip->GetResetDates();
			
			ARM_CalibMethod* firstcalibMethod = new ARM_CalibMethod(itsSwaptionCMSPF_1_2[0],ARM_ModelParamVector(),ARM_CalibMethodType::Optimize,
										 SIGMA_CSO_ARM_MAX_ITER,ARM_CalibrationTarget::PriceTarget);
			int size = resetDates->size();
			int portfolioSize = size-GetNbFixFlows();

			ARM_CalibMethod* previousCalibMethod = firstcalibMethod; 

			for (int i = 1; i<portfolioSize; i++)
			{
				ARM_CalibMethod* calibMethod_i = new ARM_CalibMethod(itsSwaptionCMSPF_1_2[i],ARM_ModelParamVector(),ARM_CalibMethodType::Optimize,
											 SIGMA_CSO_ARM_MAX_ITER,ARM_CalibrationTarget::PriceTarget);

				calibMethod_i->SetPreviousMethod(previousCalibMethod);
				previousCalibMethod = calibMethod_i;		
			} 
			
			volCalibMethod = previousCalibMethod;
		}

		return ARM_CalibMethodPtr( volCalibMethod );
}
/////////////////////////////////////////////////////////////////
///	Class  : ARM_CSOCalculator
///	Routine: computeSigmaPortfolioPrices
///	Returns: void
///	Action : create an empty calibration structure
/////////////////////////////////////////////////////////////////

void ARM_CSOCalculator::computeSigmaPortfolioPrices()
{
	ARM_BSModel* CFBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[CfModelKey]) );
	ARM_SFRM* refModel;
	refModel = dynamic_cast< ARM_SFRM* >(&*GetPricingModel());

	double tenor = ((ARM_ModelParamsSFRM*) refModel->GetModelParams())->GetIRIndex()->GetYearTerm();
	
	if (itsSigmaCalib == BootstrapCalib)
	{
		double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
		ARM_GP_Vector* resetDates = itsCalibDateStrip->GetResetDates();
		int calibSize = (*resetDates).size()-GetNbFixFlows();

		ARM_StdPortfolioPtr portfolioCMS1 = GetSigmaCalibMethod()->GetPortfolio();
		ARM_Swaption* swaptionCMS1;
		
		ARM_GP_Vector calibTimes(calibSize);
		ARM_GP_Vector volatility(calibSize);
		ARM_GP_Vector lowerVolatility(calibSize,SIGMA_LOWER_BOUND);
		ARM_GP_Vector upperVolatility(calibSize,SIGMA_UPPER_BOUND);
		for(int i =0; i<calibSize; i++)
		{
			calibTimes[i] = (*resetDates)[GetNbFixFlows()+i]-asOfDate;
			volatility[i] = CFBSModel->GetVolatility()->ComputeVolatility(calibTimes[i]/K_YEAR_LEN,tenor)/100.0;
			swaptionCMS1 = static_cast< ARM_Swaption* >(portfolioCMS1->GetAsset(i));
			swaptionCMS1->SetModel(CFBSModel);
			double swaptionCMS1_price = swaptionCMS1->ComputePrice();
			portfolioCMS1->SetPrice(swaptionCMS1_price,i);
		}
		ARM_CurveModelParam* volParam = new ARM_CurveModelParam(ARM_ModelParamType::Volatility, &volatility, &calibTimes ,"","STEPUPRIGHT",&lowerVolatility,&upperVolatility);
		
		int paramSize=GetSigmaCalibMethod()->GetCalibParams().size();

		if(paramSize == 0 )
			GetSigmaCalibMethod()->GetCalibParams().push_back(volParam);
		else
		{
			delete GetSigmaCalibMethod()->GetCalibParam(0);
			(GetSigmaCalibMethod()->GetCalibParams())[0]= volParam;
		}
	}
	else
	{
		double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
		ARM_GP_Vector* resetDates = itsCalibDateStrip->GetResetDates();
		int size = resetDates->size();
		int portfolioSize = size-GetNbFixFlows();

		double lastCalibTime = (*resetDates)[size-1]-asOfDate;
		double lastVol = CFBSModel->GetVolatility()->ComputeVolatility(lastCalibTime/K_YEAR_LEN,tenor)/100.0;
		ARM_GP_Vector breakpointTimes(1,lastCalibTime);
		ARM_GP_Vector volatility(1,lastVol);
		ARM_GP_Vector lowerVolatility(1,SIGMA_LOWER_BOUND);
		ARM_GP_Vector upperVolatility(1,SIGMA_UPPER_BOUND);

		ARM_CurveModelParam* lastVolParam = new ARM_CurveModelParam(ARM_ModelParamType::Volatility, &volatility, &breakpointTimes ,"","STEPUPRIGHT",&lowerVolatility,&upperVolatility);


		size_t sigmaIdx,paramSize=GetSigmaCalibMethod()->GetCalibParams().size();
		for(sigmaIdx=0;sigmaIdx<paramSize;++sigmaIdx)
			if(GetSigmaCalibMethod()->GetCalibParam(sigmaIdx) &&
			   (GetSigmaCalibMethod()->GetCalibParams())[sigmaIdx]->GetType() == ARM_ModelParamType::Volatility)
				break;	

		if(sigmaIdx >= paramSize || paramSize == 0)
			GetSigmaCalibMethod()->GetCalibParams().push_back(lastVolParam);
		else
		{
			delete GetSigmaCalibMethod()->GetCalibParam(sigmaIdx);
			(GetSigmaCalibMethod()->GetCalibParams())[sigmaIdx]= lastVolParam;
		}

		ARM_Swaption* swaptionCMS_k;
		ARM_StdPortfolioPtr lastPortfolioCMS_1_2 = GetSigmaCalibMethod()->GetPortfolio();
		int lastSize = lastPortfolioCMS_1_2->GetSize();
		if (itsSigmaCalib == VegaBestFitCalib)
		{
			for(int k=0;k<lastSize; k++)
			{
				swaptionCMS_k = static_cast< ARM_Swaption* >(lastPortfolioCMS_1_2->GetAsset(k));
				swaptionCMS_k->SetModel(CFBSModel);
				double swaptionCMS_k_vega = swaptionCMS_k->CptBSSensitivity(K_VEGA);
				double swaptionCMS_k_price = swaptionCMS_k->GetPrice();
				lastPortfolioCMS_1_2->SetPrice(swaptionCMS_k_price,k);
				/*
				isNotSelected = vega < VEGA_MIN_TO_SELECT * nominal && !isFreezeWeights;

				vanillaPF->SetWeight(isNotSelected ? 0.0 : weight,i);

				/// Set precision to max(1/4 vega,0.5bp)
				precision = 0.25*vega;
				minPrec = 0.00005*nominal;
				if(precision < minPrec)
					precision = minPrec;
				vanillaPF->SetPrecision(precision,i);
				*/
				lastPortfolioCMS_1_2->SetWeight(DEFAULT_WEIGHT, k);
				lastPortfolioCMS_1_2->SetPrecision(swaptionCMS_k_vega, k);
			}
		}
		else
		{
			for(int k=0;k<lastSize; k++)
			{
				swaptionCMS_k = static_cast< ARM_Swaption* >(lastPortfolioCMS_1_2->GetAsset(k));
				swaptionCMS_k->SetModel(CFBSModel);				 
				double swaptionCMS_k_price = swaptionCMS_k->ComputePrice();
				lastPortfolioCMS_1_2->SetPrice(swaptionCMS_k_price,k);
				lastPortfolioCMS_1_2->SetWeight(DEFAULT_WEIGHT, k);
				lastPortfolioCMS_1_2->SetPrecision(DEFAULT_PRECISION, k);
			}
		}

		ARM_CalibMethod* previousCalibMethod = GetSigmaCalibMethod()->GetPreviousMethod(); 

		for (int i = size-2; i>=GetNbFixFlows(); i--)
		{
			ARM_Swaption* swaptionCMS_k_i;
			double calibTime_i = (*resetDates)[i]-asOfDate;
			ARM_GP_Vector breakpointTimes_i(1,calibTime_i);
			double vol_i = CFBSModel->GetVolatility()->ComputeVolatility(calibTime_i/K_YEAR_LEN,tenor)/100.0;
			ARM_GP_Vector volatility_i(1,vol_i);
			ARM_GP_Vector lowerVolatility_i(1,SIGMA_LOWER_BOUND);
			ARM_GP_Vector upperVolatility_i(1,SIGMA_UPPER_BOUND);
			ARM_CurveModelParam* volParam_i = new ARM_CurveModelParam( ARM_ModelParamType::Volatility, &volatility_i, &breakpointTimes_i ,"","STEPUPRIGHT",&lowerVolatility_i,&upperVolatility_i);
			int paramSize_i=previousCalibMethod->GetCalibParams().size();

			if(paramSize_i == 0 )
				(previousCalibMethod->GetCalibParams()).push_back(volParam_i);
			else
			{
				delete previousCalibMethod->GetCalibParam(0);
				(previousCalibMethod->GetCalibParams())[0]= volParam_i;
			}

			ARM_StdPortfolioPtr PortfolioCMS_1_2_i = previousCalibMethod->GetPortfolio();
			int Size_i = PortfolioCMS_1_2_i->GetSize();
			if (itsSigmaCalib == VegaBestFitCalib)
			{
				for(int k=0;k<Size_i; k++)
				{
					swaptionCMS_k_i = static_cast< ARM_Swaption* >(PortfolioCMS_1_2_i->GetAsset(k));
					swaptionCMS_k_i->SetModel(CFBSModel);
					double swaptionCMS_k_vega_i = swaptionCMS_k_i->CptBSSensitivity(K_VEGA);
					double swaptionCMS_k_price_i = swaptionCMS_k_i->GetPrice();
					PortfolioCMS_1_2_i->SetPrice(swaptionCMS_k_price_i,k);
					PortfolioCMS_1_2_i->SetWeight(DEFAULT_WEIGHT, k);
					PortfolioCMS_1_2_i->SetPrecision(swaptionCMS_k_vega_i, k);
				}
			}
			else
			{
				for(int k=0;k<Size_i; k++)
				{
					swaptionCMS_k_i = static_cast< ARM_Swaption* >(PortfolioCMS_1_2_i->GetAsset(k));
					swaptionCMS_k_i->SetModel(CFBSModel);				 
					double swaptionCMS_k_price_i = swaptionCMS_k_i->ComputePrice();
					PortfolioCMS_1_2_i->SetPrice(swaptionCMS_k_price_i,k);
					PortfolioCMS_1_2_i->SetWeight(DEFAULT_WEIGHT, k);
					PortfolioCMS_1_2_i->SetPrecision(DEFAULT_PRECISION, k);
				}
			}
			previousCalibMethod = previousCalibMethod->GetPreviousMethod();
		}   
	}
}
/////////////////////////////////////////////////////////////////
///	Class  : ARM_CSOCalculator
///	Routine: CreateSOAndSWOPTPF
///	Returns: a portfolio
///	Action : create the list of diagonal swaptions
/////////////////////////////////////////////////////////////////
void ARM_CSOCalculator::CreateSOAndSWOPTPF()
{
	CreateDiagonalSwaption();
	int sizeSO = (*itsSpreadOptionPF).size();
	int sizeSWOPT = (*itsSwaptionPF).size();

	list< ARM_Security* > SOAndSwoptList;
	
    for (int i = 0; i < sizeSO; ++i)
	{
		SOAndSwoptList.push_back(static_cast<ARM_SpreadOption*>((*itsSpreadOptionPF).GetAsset(i)->Clone()));
	}

	for (i = 0; i < sizeSWOPT; ++i)
	{
		SOAndSwoptList.push_back(static_cast<ARM_Swaption*>((*itsSwaptionPF).GetAsset(i)->Clone()));
	}

	itsSOAndSWOPtPF = ARM_StdPortfolioPtr(new ARM_StdPortfolio(SOAndSwoptList));

}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CSOCalculator
///	Routine: CreateEmptyCalibration
///	Returns: void
///	Action : create an empty calibration structure
/////////////////////////////////////////////////////////////////
void ARM_CSOCalculator::CreateEmptyCalibration()
	
{
	ARM_CalibMethodPtr volCalibMethod = CreateSigmaCalibration();


	if ( (itsCalibND == ONEDCalib) || (itsMrsCalib) )
	{
		ARM_CalibMethod* ThetacalibMethod = new ARM_CalibMethod(itsSpreadOptionPF,ARM_ModelParamVector(),ARM_CalibMethodType::Optimize,
											 CORREL_CSO_ARM_MAX_ITER,ARM_CalibrationTarget::PriceTarget);

		if(itsMrsCalib)
		{
			CreateDiagonalSwaption();
			ARM_CalibMethod* MrscalibMethod = new ARM_CalibMethod(itsSwaptionPF,ARM_ModelParamVector(),ARM_CalibMethodType::Optimize,
											 MRS_CSO_ARM_MAX_ITER,ARM_CalibrationTarget::PriceTarget);

			MrscalibMethod->SetlinkedMethod((ARM_CalibMethod*)(&(*volCalibMethod))->Clone());
			ThetacalibMethod->SetlinkedMethod(MrscalibMethod);
		}
		else
		{
			ThetacalibMethod->SetlinkedMethod((ARM_CalibMethod*)(&(*volCalibMethod))->Clone());
		}

		SetCalibMethod(ARM_CalibMethodPtr(ThetacalibMethod));

	}
	else
	{
		CreateSOAndSWOPTPF();

		ARM_CalibMethod* ThetaMrscalibMethod = new ARM_CalibMethod(itsSOAndSWOPtPF,ARM_ModelParamVector(),ARM_CalibMethodType::Optimize,
											 CORREL_CSO_ARM_MAX_ITER,ARM_CalibrationTarget::PriceTarget);

		ThetaMrscalibMethod->SetlinkedMethod((ARM_CalibMethod*)(&(*volCalibMethod))->Clone());

		SetCalibMethod(ARM_CalibMethodPtr(ThetaMrscalibMethod));
	}
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CSOCalculator
///	Routine: CreateAndSetCalibration
///	Returns: void
///	Action : create the calibration
/////////////////////////////////////////////////////////////////
void ARM_CSOCalculator::CreateAndSetCalibration()
{
	CreateEmptyCalibration();

	computeSigmaPortfolioPrices();

	if ( (itsCalibND == ONEDCalib) || (itsMrsCalib) )
	{
		computeSOPortfolioPrices();

		if(itsMrsCalib)
		{
			computeDIAGSwaptionPortfolioPrices();
		}
	}
	else
	{
		//initial model
		computeSOAndSWOPTPortfolioPrices();
	}
	
}

////////////////////////////////////////////////////////////////
///	Class  : ARM_CSOCalculator
///	Routine: computeDIAGSwaptionPortfolioPrices
///	Returns: nothing
///	Action : compute market target prices of the swaption portfolio
/////////////////////////////////////////////////////////////////
void ARM_CSOCalculator::computeSOAndSWOPTPortfolioPrices()
{
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	
	ARM_BSModel* SOBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[SoModelKey]) );
	ARM_BSModel* OswBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswModelKey]) );
	//2) pricer le portfolio SOAndSWOPT et Calculer VEga et Cega
	ARM_StdPortfolioPtr soAndSwoptPF = GetCalibMethod()->GetPortfolio();
	int sizeSO = (*itsSpreadOptionPF).size();
	int sizeSWOPT = (*itsSwaptionPF).size();
	ARM_SpreadOption* so;
	ARM_Swaption* swopt;
    for (int k = 0; k < sizeSO; ++k)
	{
		so = static_cast< ARM_SpreadOption* >(soAndSwoptPF->GetAsset(k));
		so->SetModel(SOBSModel);
		double SO_k_vega = so->ComputeSensitivity(K_VEGA);
		double SO_k_correl = so->ComputeSensitivity(K_CORREL);
		double weight = fabs(SO_k_vega+SO_k_correl);
		double SO_k_price = so->ComputePrice();
		soAndSwoptPF->SetPrice(SO_k_price,k);
		soAndSwoptPF->SetWeight(itsModelDatas[4], k);
		soAndSwoptPF->SetPrecision(weight, k);
		itsSpreadOptionPF->SetPrice(SO_k_price,k);
		itsSpreadOptionPF->SetWeight(DEFAULT_WEIGHT, k);
		itsSpreadOptionPF->SetPrecision(weight, k);
	}
	for (k = sizeSO; k < (sizeSWOPT+sizeSO); ++k)
	{
		swopt = static_cast< ARM_Swaption* >(soAndSwoptPF->GetAsset(k));
		swopt->SetModel(OswBSModel);
		double swopt_k_vega = swopt->ComputeSensitivity(K_VEGA);
		double swopt_k_price = swopt->ComputePrice();
		soAndSwoptPF->SetPrice(swopt_k_price,k);
		soAndSwoptPF->SetWeight(DEFAULT_WEIGHT, k);
		soAndSwoptPF->SetPrecision(swopt_k_vega, k);
	}
	if(itsUpdateInitCalib)
	{
		//1) lancer une calib theta et sigma
		ARM_CalibMethod* initSigmaCalibmethod = (ARM_CalibMethod*) GetSigmaCalibMethod()->Clone();
		ARM_Date cmsDate1(GetEndDate());
		if (IsCMSIndex((ARM_INDEX_TYPE)GetCMSLong()))
			cmsDate1.AddMonths( (GetCMSLong()-K_CMS1+1)*12 );

		ARM_DateStrip datestrip(GetStartDate(),cmsDate1,CORREL_FREQ,GetCpnDaycount(),
							GETDEFAULTVALUESTR,
							K_MOD_FOLLOWING,K_ADJUSTED,K_SHORTSTART,
							GETDEFAULTVALUE,CORREL_FREQ);

		ARM_GP_Vector lowerBound(1,CORREL_LOWER_BOUND);
		ARM_GP_Vector upperBound(1,CORREL_UPPER_BOUND);

		ARM_CalibMethod* ThetacalibMethod = new ARM_CalibMethod(itsSpreadOptionPF,ARM_ModelParamVector(),ARM_CalibMethodType::Optimize,
												 CORREL_CSO_ARM_MAX_ITER,ARM_CalibrationTarget::PriceTarget);

		ARM_CorrelTrigoMatParam* correlparam = new ARM_CorrelTrigoMatParam(CORREL_DEFAULT_INITVALUE,
																datestrip,
																asOfDate,
																"LINEAR",
																&lowerBound,
																&upperBound);

		ThetacalibMethod->GetCalibParams().push_back(correlparam);
			
		//pricer le SO

		ThetacalibMethod->SetlinkedMethod(initSigmaCalibmethod);

		ARM_SFRM* refModel = dynamic_cast< ARM_SFRM* >(GetPricingModel()->Clone());

		//initSFRM = calibrate(SFRM,ThetacalibMethod)
		ThetacalibMethod->Calibrate(refModel);
				
		itsTWODInitSigmaParam = ARM_ModelParamPtr(static_cast<ARM_ModelParam*>(refModel->GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).Clone()));
		itsTWODInitThetaParam = ARM_ModelParamPtr(static_cast<ARM_ModelParam*>(refModel->GetModelParams()->GetModelParam(ARM_ModelParamType::Correlation).Clone()));
	
		delete refModel;
		delete ThetacalibMethod;
		//initSFRM getsigmaParam, getthetaparam
		//itsTWODInitSigmaParam, itsTWODInitThetaParam;		
		itsUpdateInitCalib = false;
	}

	// model param theta et le sigma et Mrs
	ARM_GP_Vector initTimes(1,0.0);
	ARM_GP_Vector initMRS(1,MRS_DEFAULT_VALUE);
	ARM_GP_Vector mrsLowerBound(1,MRS_LOWER_BOUND);
	ARM_GP_Vector mrsUpperBound(1,MRS_UPPER_BOUND);

	ARM_CurveModelParam* mrs = new ARM_CurveModelParam(
		ARM_ModelParamType::MeanReversion,
		&initMRS,
		&initTimes,
		"MEANREVERSION",
		"STEPUPRIGHT",
		&mrsLowerBound,
		&mrsUpperBound);


	ARM_CalibMethod* twoDimCalibmethod = &*GetCalibMethod();
	int paramSize=twoDimCalibmethod->GetCalibParams().size();

	if(paramSize == 0 )//CLEAN
	{
		(twoDimCalibmethod->GetCalibParams()).push_back((ARM_CurveModelParam*)(((ARM_CorrelTrigoMatParam*)(&*itsTWODInitThetaParam))->Clone()));
		(twoDimCalibmethod->GetCalibParams()).push_back(mrs);
	}
	else //CLEAN
	{
		delete twoDimCalibmethod->GetCalibParam(0);
		delete twoDimCalibmethod->GetCalibParam(1);		
		(twoDimCalibmethod->GetCalibParams())[0]= (ARM_CurveModelParam*)(((ARM_CorrelTrigoMatParam*)(&*itsTWODInitThetaParam))->Clone());
		(twoDimCalibmethod->GetCalibParams())[1]= mrs;
		//mrs 
	}


	ARM_CalibMethod* sigmaCalibmethod  = GetSigmaCalibMethod();
	
	ARM_GP_Vector* resetDates = itsCalibDateStrip->GetResetDates();
	int size = resetDates->size();
	int portfolioSize = size-GetNbFixFlows();
	ARM_Curve* sigmaCurve = ((ARM_CurveModelParam*) (&*itsTWODInitSigmaParam))->GetCurve();
	
	double lastCalibTime = (*resetDates)[size-1]-asOfDate;
	double lastVol = sigmaCurve->Interpolate(lastCalibTime);
	ARM_GP_Vector breakpointTimes(1,lastCalibTime);
	ARM_GP_Vector volatility(1,lastVol);
	ARM_GP_Vector lowerVolatility(1,SIGMA_LOWER_BOUND);
	ARM_GP_Vector upperVolatility(1,SIGMA_UPPER_BOUND);	
	ARM_CurveModelParam* lastVolParam = new ARM_CurveModelParam(ARM_ModelParamType::Volatility, &volatility, &breakpointTimes ,"","STEPUPRIGHT",&lowerVolatility,&upperVolatility);

	paramSize=GetSigmaCalibMethod()->GetCalibParams().size();

	if(paramSize == 0 )
		GetSigmaCalibMethod()->GetCalibParams().push_back(lastVolParam);
	else
	{
		delete GetSigmaCalibMethod()->GetCalibParam(0);
		(GetSigmaCalibMethod()->GetCalibParams())[0]= lastVolParam;
	}

	ARM_CalibMethod* previousCalibMethod = GetSigmaCalibMethod()->GetPreviousMethod(); 

	for (int i = size-2; i>=GetNbFixFlows(); i--)
	{
		double calibTime_i = (*resetDates)[i]-asOfDate;
		ARM_GP_Vector breakpointTimes_i(1,calibTime_i);
		double vol_i = sigmaCurve->Interpolate(calibTime_i);
		ARM_GP_Vector volatility_i(1,vol_i);
		ARM_GP_Vector lowerVolatility_i(1,SIGMA_LOWER_BOUND);
		ARM_GP_Vector upperVolatility_i(1,SIGMA_UPPER_BOUND);
		ARM_CurveModelParam* volParam_i = new ARM_CurveModelParam( ARM_ModelParamType::Volatility, &volatility_i, &breakpointTimes_i ,"","STEPUPRIGHT",&lowerVolatility_i,&upperVolatility_i);
		int paramSize_i=previousCalibMethod->GetCalibParams().size();

		if(paramSize_i == 0 )
			(previousCalibMethod->GetCalibParams()).push_back(volParam_i);
		else
		{
			delete previousCalibMethod->GetCalibParam(0);
			(previousCalibMethod->GetCalibParams())[0]= volParam_i;
		}

		previousCalibMethod = previousCalibMethod->GetPreviousMethod();
	}	
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CSOCalculator
///	Routine: GetOSWCalibMethod
///	Returns: ARM_CalibMethod
///	Action : get calibration method for diagonal swaptions
/////////////////////////////////////////////////////////////////
ARM_CalibMethod* ARM_CSOCalculator::GetSWOPTCalibMethod() const
{
	if(itsMrsCalib == false)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : No method available because calibration is off");

    return GetCalibMethod()->GetlinkedMethod();
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CSOCalculator
///	Routine: GetSOPortfolio
///	Returns: ARM_Portfolio
///	Action : get the spreadoption calibration portfolio
/////////////////////////////////////////////////////////////////
const ARM_StdPortfolioPtr ARM_CSOCalculator::GetSOPortfolio() const
{
	if( GetCalibMethod() == ARM_CalibMethodPtr(NULL) ||
		GetCalibMethod()->GetPortfolio() == ARM_StdPortfolioPtr(NULL) )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Optimisation calib method or portfolio not found");

		return GetCalibMethod()->GetPortfolio();
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CSOCalculator
///	Routine: GetOSWPortfolio
///	Returns: ARM_Portfolio
///	Action : get the diagonal swaptions calibration portfolio
/////////////////////////////////////////////////////////////////
const ARM_StdPortfolioPtr ARM_CSOCalculator::GetOSWPortfolio() const
{
	if( (itsCalibND == TWODCalib) || (itsMrsCalib == false) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : No method available because calibration is off");

	return GetCalibMethod()->GetlinkedMethod()->GetPortfolio();
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CSOCalculator
///	Routine: GetOSWPortfolio
///	Returns: ARM_Portfolio
///	Action : get the diagonal swaptions calibration portfolio
/////////////////////////////////////////////////////////////////
const ARM_StdPortfolioPtr ARM_CSOCalculator::GetCMSLONGPortfolio() const
{
	list< ARM_Security* > swaptionList;

	if (itsSigmaCalib == BootstrapCalib)
		return GetSigmaCalibMethod()->GetPortfolio();
	ARM_CalibMethod* sigmaCalibMethod = GetSigmaCalibMethod();
	while (sigmaCalibMethod)
	{
		ARM_StdPortfolioPtr tmpPf = sigmaCalibMethod->GetPortfolio();
		
		swaptionList.push_back(static_cast<ARM_Swaption*>((*tmpPf).GetAsset(0)->Clone()));
		sigmaCalibMethod = sigmaCalibMethod->GetPreviousMethod();
	}

	return ARM_StdPortfolioPtr(new ARM_StdPortfolio(swaptionList));
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CSOCalculator
///	Routine: GetOSWPortfolio
///	Returns: ARM_Portfolio
///	Action : get the diagonal swaptions calibration portfolio
/////////////////////////////////////////////////////////////////
const ARM_StdPortfolioPtr ARM_CSOCalculator::GetCMSSHORTPortfolio() const
{
	list< ARM_Security* > swaptionList;

	if (itsSigmaCalib == BootstrapCalib)
		return GetSigmaCalibMethod()->GetPortfolio();
	ARM_CalibMethod* sigmaCalibMethod = GetSigmaCalibMethod();
	while (sigmaCalibMethod)
	{
		ARM_StdPortfolioPtr tmpPf = sigmaCalibMethod->GetPortfolio();
		
		swaptionList.push_back(static_cast<ARM_Swaption*>((*tmpPf).GetAsset(1)->Clone()));
		sigmaCalibMethod = sigmaCalibMethod->GetPreviousMethod();
	}

	return ARM_StdPortfolioPtr(new ARM_StdPortfolio(swaptionList));
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CSOCalculator
///	Routine: UpdateModel
///	Returns: void
///	Action : update the model datas w.r.t. MarketDataManager
/////////////////////////////////////////////////////////////////
void ARM_CSOCalculator::UpdateModel()
{
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CSOCalculator
///	Routine: UpdateCalibration
///	Returns: void
///	Action : update the calibration datas w.r.t. MarketDataManager
///          The context is an hedge ratios computation
/////////////////////////////////////////////////////////////////
void ARM_CSOCalculator::UpdateCalibration(bool isUpdateStrike)
{
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CSOCalculator
///	Routine: Calibrate
///	Returns: void
///	Action : Calibrate the model of the calculator.
/////////////////////////////////////////////////////////////////
void ARM_CSOCalculator::Calibrate()
{
	/// Volatility bootstrapping with possible MRS optimisation
    /// FIX FIX : ugly test waiting for the generic calibration extended for hybrid models
    /// with an equivalent modelMap on modelParam
	GetCalibMethod()->Calibrate(&*GetPricingModel());   

}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CSOCalculator
///	Routine: Price
///	Returns: a double
///	Action : price the TARN deal.
/////////////////////////////////////////////////////////////////
double ARM_CSOCalculator::Price()
{
	CalibrateAndTimeIt();

    /// Price the implicit product according to internal flag
	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();

	CC_Ostringstream evalDateDesc;
    evalDateDesc << CC_NS(std,fixed) << asOfDate.GetJulian();
    string evalDateStr(evalDateDesc.str());

	ARM_ZeroCurve* zeroCurve = (ARM_ZeroCurve*) GetMktDataManager()->GetData(GetKeys()[YcKey]);

    for (size_t i = 0; i < NbProductsToPrice; ++i)
	{
		if (itsProductsToPrice[i])
		{
			/// Select the right column that describes the product
			CSOColAlias prodIdx;

			if(i == CSOPrice)
				prodIdx = CSO;							
			else if(i == structPrice)
				prodIdx = SpreadOptionLet;
			
			ARM_GenSecurityPtr genSec = GetGenSecurity();
			ARM_DealDescription dealDesc = GetGenSecurity()->GetDealDescription();

			/// Build the sub-generic security
			size_t startRowIdx = 1; /// start with the first actual row (=> skip column names)
			size_t endRowIdx = dealDesc.GetRowsNb() - 1;
			ARM_DealDescriptionPtr subDealDesc = dealDesc.GetSubDescription(startRowIdx,endRowIdx,prodIdx+1);
			ARM_GenSecurityPtr subGenSec = ARM_GenSecurityPtr(new ARM_GenSecurity(subDealDesc,genSec->GetPayModelName()));

			//genSec = ARM_GenSecurityPtr(new ARM_GenSecurity(subDealDesc,genSec->GetPayModelName()));

			ARM_GenPricer* genPricer = new ARM_GenPricer( &*subGenSec,&*GetPricingModel());

			ARM_AutoCleaner<ARM_GenPricer> HoldGP(genPricer );

			double price = genPricer->Price();
			
			if(i == CSOPrice)
			{
				itsCSOPrice = price;
			}
			else if(i == structPrice)
			{
				itsStructPrice = price;
			}
		}
	}

	itsHasBeenPriced = true;
    
    return itsCSOPrice;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CSOCalculator
///	Routine: CheckData & CheckMktData
///	Returns: void
///	Action : check if TARN datas are consistent
/////////////////////////////////////////////////////////////////
void ARM_CSOCalculator::CheckData()
{
}

void ARM_CSOCalculator::CheckMktData()
{
}



CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
