/*!
 *
 * Copyright (c) IXIS CIB August 2006 Paris
 *
 *	\file tarnfxcalculator.cpp
 *
 *  \brief file for the TARN FX Calculator
 *	\author  A. Lekrafi
 *	\version 1.0
 *	\date August 2006
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpcalculators/tarncalculatorindian.h"
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

const string ARM_TARNCalculatorIndian::TARNIndianColNamesTable[] =
{
	"ResetDate",
	"PayDate",
	"BarrierUp",
	"BarrierDown",
	"Strike",
	"Notional",
	"Target",
	"Fees",
	"SpotFwd",
	"TwiceCallPut",
	"CallMinusCall",
	"PutMinusPut",
	"Epsilon",
	"CallDigital",
	"PutDigital",
	"IsAlive",
	"Coupon",
	"SumCoupon",
	"PaidCoupon",
	"RealCoupon",
	"Proba",
	"TARN"
};


const int ARM_TARNCalculatorIndian::ProductToPriceColumns[] =
{
	ARM_TARNCalculatorIndian::SpotFwd,
	ARM_TARNCalculatorIndian::TwiceCallPut,
	ARM_TARNCalculatorIndian::CallMinusCall,
	ARM_TARNCalculatorIndian::PutMinusPut,
	ARM_TARNCalculatorIndian::CallDigital,
	ARM_TARNCalculatorIndian::PutDigital,
	ARM_TARNCalculatorIndian::IsAlive,
	ARM_TARNCalculatorIndian::Coupon,
	ARM_TARNCalculatorIndian::PaidCoupon,
	ARM_TARNCalculatorIndian::RealCoupon,
	ARM_TARNCalculatorIndian::Proba,
	ARM_TARNCalculatorIndian::TARN,
};


/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculatorIndian
///	Routine: Copy constructor
///	Returns: a built object
///	Action : builds the object
/////////////////////////////////////////////////////////////////
ARM_TARNCalculatorIndian::ARM_TARNCalculatorIndian( const ARM_TARNCalculatorIndian& rhs )
:ARM_TARNFXCalculator(rhs)
{
	itsBarrierUp	= rhs.itsBarrierUp;
	itsBarrierDown	= rhs.itsBarrierDown;
	itsEpsilon		= rhs.itsEpsilon;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculatorIndian
///	Routine: Destructor
///	Returns: void
///	Action : destroys the object
/////////////////////////////////////////////////////////////////
ARM_TARNCalculatorIndian::~ARM_TARNCalculatorIndian()
{
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculatorIndian
///	Routine: Constructor
///	Returns: void
///	Action : builds the object (stand alone version)
/////////////////////////////////////////////////////////////////
ARM_TARNCalculatorIndian::ARM_TARNCalculatorIndian ( const ARM_Date& asOfDate,
													 const ARM_Date& startDate,
													 const ARM_Date& endDate,
													 const ARM_Currency& DomCcy,
													 const vector<ARM_Currency>& ForCcy,
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
													 const ARM_Curve& barrierUp,
													 const ARM_Curve& barrierDown,
													 double target,
													 double epsilon,
													 const ARM_Curve& fees,
													 int intermediatePrices,
													 const ARM_IntVector& productsToPrice)
:
	ARM_TARNFXCalculator(asOfDate,
						 startDate,
						 endDate,
						 DomCcy,
						 ForCcy,
						 DomCcy,
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
						 vector<ARM_Curve>(),
						 vector<ARM_Curve>(),
						 ARM_Curve(),
						 ARM_Curve(),
						 vector<ARM_Curve>(),
						 K_DEF_FREQ,
						 KNOBASE,
						 ARM_Curve(),
						 ARM_Curve(),
						 target,
						 "",
						 0,
						 0,
						 ARM_GP_Vector(),
						 fees,
						 intermediatePrices,
						 productsToPrice),
	itsBarrierUp(barrierUp),
	itsBarrierDown(barrierDown),
	itsEpsilon(epsilon)
{
	string payModelName = YC_BASIS_KEY_NAME + string(itsDomesticCcy.GetCcyName());

	itsColumnsToPrice = ProductsToPriceColumnNames();

    /// Create the Generic Security
	CreateAndSetDealDescription(payModelName, itsColumnsToPrice, CreateCstManager(), false);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculatorIndian
///	Routine: CheckData
///	Returns: void
///	Action : check if TARN FX data are consistent
/////////////////////////////////////////////////////////////////
void ARM_TARNCalculatorIndian::CheckData()
{
	if (itsStartDate >= itsEndDate)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : The start date frequency should be before the end date.");
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculatorIndian
///	Routine: CheckMktData
///	Returns: void
///	Action : check if TARN FX market data are consistent
/////////////////////////////////////////////////////////////////
void ARM_TARNCalculatorIndian::CheckMktData()
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
/*
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
*/
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

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculatorIndian
///	Routine: ColumnNames
///	Returns: ARM_RowInfo
///	Action : create the column names of the deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_TARNCalculatorIndian::ColumnNames() const
{
    size_t	colNamesSize = sizeof(TARNIndianColNamesTable)/sizeof(TARNIndianColNamesTable[0]);
    vector<string>				colNamesVec(colNamesSize);
    vector<ARM_GP_VALUE_TYPE>	colTypeVec(colNamesSize, ARM_STRING); 

    for(size_t i=0; i<colNamesSize; ++i)
        colNamesVec[i] = TARNIndianColNamesTable[i];

    ARM_RowInfo	rowInfo(colNamesVec, colTypeVec);

    return rowInfo;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculatorIndian
///	Routine: ProductsToPriceColumnNames
///	Returns: ARM_StringVector
///	Action : Create Products to price column names
/////////////////////////////////////////////////////////////////
ARM_StringVector ARM_TARNCalculatorIndian::ProductsToPriceColumnNames()
{
	ARM_StringVector columnNames;
	columnNames.reserve(NbProductsToPrice);

    for(size_t i = 0; i < NbProductsToPrice; ++i)
		if ((itsProductsToPrice.size() > i) && itsProductsToPrice[i])
			columnNames.push_back( TARNIndianColNamesTable[ ARM_TARNCalculatorIndian::ProductToPriceColumns[i] ] );

	return columnNames;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculatorIndian
///	Routine: InitPriceableColumns
///	Returns: void
///	Action : initialise to 0 column to be able to sum
/////////////////////////////////////////////////////////////////
void ARM_TARNCalculatorIndian::InitPriceableColumns(vector<string>& rowDescVec, vector<ARM_GP_VALUE_TYPE>& rowTypeVec) const
{
    string	zeroValue("0");

	rowDescVec[ARM_TARNCalculatorIndian::SpotFwd] = zeroValue;
    rowTypeVec[ARM_TARNCalculatorIndian::SpotFwd] = ARM_DOUBLE;

	rowDescVec[ARM_TARNCalculatorIndian::TwiceCallPut] = zeroValue;
    rowTypeVec[ARM_TARNCalculatorIndian::TwiceCallPut] = ARM_DOUBLE;

	rowDescVec[ARM_TARNCalculatorIndian::CallMinusCall] = zeroValue;
    rowTypeVec[ARM_TARNCalculatorIndian::CallMinusCall] = ARM_DOUBLE;

	rowDescVec[ARM_TARNCalculatorIndian::PutMinusPut] = zeroValue;
    rowTypeVec[ARM_TARNCalculatorIndian::PutMinusPut] = ARM_DOUBLE;

	rowDescVec[ARM_TARNCalculatorIndian::CallDigital] = zeroValue;
    rowTypeVec[ARM_TARNCalculatorIndian::CallDigital] = ARM_DOUBLE;

	rowDescVec[ARM_TARNCalculatorIndian::PutDigital] = zeroValue;
    rowTypeVec[ARM_TARNCalculatorIndian::PutDigital] = ARM_DOUBLE;

	rowDescVec[ARM_TARNCalculatorIndian::IsAlive] = zeroValue;
    rowTypeVec[ARM_TARNCalculatorIndian::IsAlive] = ARM_DOUBLE;

	rowDescVec[ARM_TARNCalculatorIndian::Coupon] = zeroValue;
    rowTypeVec[ARM_TARNCalculatorIndian::Coupon] = ARM_DOUBLE;

	rowDescVec[ARM_TARNCalculatorIndian::PaidCoupon] = zeroValue;
    rowTypeVec[ARM_TARNCalculatorIndian::PaidCoupon] = ARM_DOUBLE;

	rowDescVec[ARM_TARNCalculatorIndian::RealCoupon] = zeroValue;
    rowTypeVec[ARM_TARNCalculatorIndian::RealCoupon] = ARM_DOUBLE;

	rowDescVec[ARM_TARNCalculatorIndian::Proba] = zeroValue;
    rowTypeVec[ARM_TARNCalculatorIndian::Proba] = ARM_DOUBLE;

	rowDescVec[ARM_TARNCalculatorIndian::TARN] = zeroValue;
    rowTypeVec[ARM_TARNCalculatorIndian::TARN] = ARM_DOUBLE;
}


//////////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculatorIndian
///	Routine: DatesStructure
///	Returns: ARM_DateStripCombiner
///	Action : create the list of all event/notice dates of the TARN FX.
//////////////////////////////////////////////////////////////////////
ARM_DateStripCombiner ARM_TARNCalculatorIndian::DatesStructure() const
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
	}

    /// Merge schedules on "ResetDate"
    ARM_DateStripVector SchedVect(1, NULL);
    SchedVect[RF_CPNFX_SCHED] = &*itsCpnDateStrip;

    ARM_DateStripCombiner EventSchedule(SchedVect, "ResetDate");

	// Memorize first event index
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	ARM_VectorPtr eventDates = EventSchedule.GetMergeData();

	int	firstEventIdx = 0;
	while( (firstEventIdx < eventDates->size()) && ((*eventDates)[firstEventIdx] < asOfDate) )
		++firstEventIdx;

	const_cast<ARM_TARNCalculatorIndian*>(this)->itsFirstEventIdx = firstEventIdx;

    return	EventSchedule;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculatorIndian
///	Routine: MiddleRows
///	Returns: ARM_RowInfo
///	Action : create a row of a deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_TARNCalculatorIndian::MiddleRows( size_t eventIdx, const ARM_DateStripCombiner& datesStructure ) const
{
	ARM_DateStripPtr	vDateStrip = datesStructure.GetDateStrip(0);
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	size_t	eventSize = vDateStrip->GetResetDates()->size();
	size_t	descSize = sizeof(TARNIndianColNamesTable)/sizeof(TARNIndianColNamesTable[0]);

	vector<string>				rowDescVec(descSize);
	vector<ARM_GP_VALUE_TYPE>	rowTypeVec(descSize, ARM_MISSING_TYPE); 

    /// Set default 0 value for each column to be able to sum it
    InitPriceableColumns(rowDescVec, rowTypeVec);
	int i = 0;

	bool isFirstCpnEvent = ((itsCpnTiming == K_ARREARS)&&(itsFirstEventIdx==0)?(eventIdx == 1):(eventIdx == itsFirstEventIdx));

    CC_Ostringstream	resetDateDesc;
	double vResetDate = (*(vDateStrip->GetResetDates()))[eventIdx];
    resetDateDesc << CC_NS(std, fixed) << vResetDate;
    rowDescVec[TARNIndianColAlias::ResetDate] = resetDateDesc.str();
    rowTypeVec[TARNIndianColAlias::ResetDate] = ARM_DATE_TYPE;

	CC_Ostringstream	payDateDesc;
	double vPayDate = (*(vDateStrip->GetPaymentDates()))[eventIdx];
	payDateDesc << CC_NS(std, fixed) << vPayDate;
	rowDescVec[TARNIndianColAlias::PayDate] = payDateDesc.str();
	rowTypeVec[TARNIndianColAlias::PayDate] = ARM_DATE_TYPE;

	CC_Ostringstream	barrierUpDesc;
	barrierUpDesc << CC_NS(std, fixed) << itsBarrierUp.Interpolate(vResetDate-asOfDate);
	rowDescVec[TARNIndianColAlias::BarrierUp] = barrierUpDesc.str();
	rowTypeVec[TARNIndianColAlias::BarrierUp] = ARM_DOUBLE;

	CC_Ostringstream	barrierDownDesc;
	barrierDownDesc << CC_NS(std, fixed) << itsBarrierDown.Interpolate(vResetDate-asOfDate);
	rowDescVec[TARNIndianColAlias::BarrierDown] = barrierDownDesc.str();
	rowTypeVec[TARNIndianColAlias::BarrierDown] = ARM_DOUBLE;

	CC_Ostringstream	targetRedemptionDesc; 
	targetRedemptionDesc << CC_NS(std, fixed) << itsTargetRedemption;
	rowDescVec[TARNIndianColAlias::Target] = targetRedemptionDesc.str();
	rowTypeVec[TARNIndianColAlias::Target] = ARM_DOUBLE;

	CC_Ostringstream	notionalDesc;
	notionalDesc << CC_NS(std, fixed) << itsCpnNominalCv.Interpolate(vPayDate-asOfDate);
	rowDescVec[TARNIndianColAlias::Notional] = notionalDesc.str();
	rowTypeVec[TARNIndianColAlias::Notional] = ARM_DOUBLE;

	CC_Ostringstream	feesDesc;
	feesDesc << CC_NS(std, fixed) << itsFees.Interpolate(vPayDate-asOfDate);
	rowDescVec[TARNIndianColAlias::Fees] = feesDesc.str();
	rowTypeVec[TARNIndianColAlias::Fees] = ARM_DOUBLE;

	CC_Ostringstream	spotFwdDesc;
	spotFwdDesc << CC_NS(std, fixed) << "SPOT(FOREX_" << itsForeignCcy[0].GetCcyName() << "/" << itsDomesticCcy.GetCcyName() << ")";
	rowDescVec[TARNIndianColAlias::SpotFwd] = spotFwdDesc.str();
	rowTypeVec[TARNIndianColAlias::SpotFwd] = ARM_DOUBLE;

	CC_Ostringstream	twiceCallPutDesc;
	twiceCallPutDesc	<< "2*(Max(" << TARNIndianColNamesTable[TARNIndianColAlias::SpotFwd] << "[i]"
						<< "-" << TARNIndianColNamesTable[TARNIndianColAlias::BarrierUp] << "[i],0.0)"
						<< "+Max(" << TARNIndianColNamesTable[TARNIndianColAlias::BarrierDown] << "[i]"
						<< "-" << TARNIndianColNamesTable[TARNIndianColAlias::SpotFwd] << "[i],0.0))";
	rowDescVec[TARNIndianColAlias::TwiceCallPut] = twiceCallPutDesc.str();
	rowTypeVec[TARNIndianColAlias::TwiceCallPut] = ARM_STRING;

	CC_Ostringstream	callMinusCallDesc;
	callMinusCallDesc	<< "Max(" << TARNIndianColNamesTable[TARNIndianColAlias::SpotFwd] << "[i]"
						<< "-" << TARNIndianColNamesTable[TARNIndianColAlias::BarrierDown] << "[i],0.0)"
						<< "-Max(" << TARNIndianColNamesTable[TARNIndianColAlias::SpotFwd] << "[i]"
						<< "-" << TARNIndianColNamesTable[TARNIndianColAlias::Strike] << "[i],0.0)";
	rowDescVec[TARNIndianColAlias::CallMinusCall] = callMinusCallDesc.str();
	rowTypeVec[TARNIndianColAlias::CallMinusCall] = ARM_STRING;

	CC_Ostringstream	putMinusPutDesc;
	putMinusPutDesc	<< "Max(" << TARNIndianColNamesTable[TARNIndianColAlias::BarrierUp] << "[i]"
					<< "-" << TARNIndianColNamesTable[TARNIndianColAlias::SpotFwd] << "[i],0.0)"
					<< "-Max(" << TARNIndianColNamesTable[TARNIndianColAlias::Strike] << "[i]"
					<< "-" << TARNIndianColNamesTable[TARNIndianColAlias::SpotFwd] << "[i],0.0)";
	rowDescVec[TARNIndianColAlias::PutMinusPut] = putMinusPutDesc.str();
	rowTypeVec[TARNIndianColAlias::PutMinusPut] = ARM_STRING;

	CC_Ostringstream	epsilonDesc;
	epsilonDesc	<< "0.01*" << TARNIndianColNamesTable[TARNIndianColAlias::Strike] << "[i]";
	rowDescVec[TARNIndianColAlias::Epsilon] = epsilonDesc.str();
	rowTypeVec[TARNIndianColAlias::Epsilon] = ARM_STRING;

	CC_Ostringstream	callDigitalDesc;
	callDigitalDesc	<< "(" << TARNIndianColNamesTable[TARNIndianColAlias::Strike] << "[i]"
					<< "-" << TARNIndianColNamesTable[TARNIndianColAlias::BarrierDown] << "[i])"
					<< "*(Max(" << TARNIndianColNamesTable[TARNIndianColAlias::SpotFwd] << "[i]"
					<< "-(" << TARNIndianColNamesTable[TARNIndianColAlias::Strike] << "[i]"
					<< "-" << TARNIndianColNamesTable[TARNIndianColAlias::Epsilon] << "[i]),0.0)"
					<< "-Max(" << TARNIndianColNamesTable[TARNIndianColAlias::SpotFwd] << "[i]"
					<< "-(" << TARNIndianColNamesTable[TARNIndianColAlias::Strike] << "[i]"
					<< "+" << TARNIndianColNamesTable[TARNIndianColAlias::Epsilon] << "[i]),0.0))"
					<< "/(2*" << TARNIndianColNamesTable[TARNIndianColAlias::Epsilon] << "[i])";
	rowDescVec[TARNIndianColAlias::CallDigital] = callDigitalDesc.str();
	rowTypeVec[TARNIndianColAlias::CallDigital] = ARM_STRING;

	CC_Ostringstream	putDigitalDesc;
	putDigitalDesc	<< "(" << TARNIndianColNamesTable[TARNIndianColAlias::BarrierUp] << "[i]"
					<< "-" << TARNIndianColNamesTable[TARNIndianColAlias::Strike] << "[i])"
					<< "*(Max((" << TARNIndianColNamesTable[TARNIndianColAlias::Strike] << "[i]"
					<< "+" << TARNIndianColNamesTable[TARNIndianColAlias::Epsilon] << "[i])"
					<< "-" << TARNIndianColNamesTable[TARNIndianColAlias::SpotFwd] << "[i],0.0)"
					<< "-Max((" << TARNIndianColNamesTable[TARNIndianColAlias::Strike] << "[i]"
					<< "-" << TARNIndianColNamesTable[TARNIndianColAlias::Epsilon] << "[i])"
					<< "-" << TARNIndianColNamesTable[TARNIndianColAlias::SpotFwd] << "[i],0.0))"
					<< "/(2*" << TARNIndianColNamesTable[TARNIndianColAlias::Epsilon] << "[i])";
	rowDescVec[TARNIndianColAlias::PutDigital] = putDigitalDesc.str();
	rowTypeVec[TARNIndianColAlias::PutDigital] = ARM_STRING;

	CC_Ostringstream	couponDesc;
	couponDesc	<< TARNIndianColNamesTable[TARNIndianColAlias::CallMinusCall] << "[i]"
				<< "-" << TARNIndianColNamesTable[TARNIndianColAlias::CallDigital] << "[i]"
				<< "+" << TARNIndianColNamesTable[TARNIndianColAlias::PutMinusPut] << "[i]"
				<< "-" << TARNIndianColNamesTable[TARNIndianColAlias::PutDigital] << "[i]";
	rowDescVec[TARNIndianColAlias::Coupon] = couponDesc.str();
	rowTypeVec[TARNIndianColAlias::Coupon] = ARM_STRING;

	CC_Ostringstream	sumCouponDesc;
	if (isFirstCpnEvent)
		sumCouponDesc	<< TARNIndianColNamesTable[TARNIndianColAlias::Coupon] << "[i]";
	else
		sumCouponDesc	<< TARNIndianColNamesTable[TARNIndianColAlias::SumCoupon] << "[i-1]"
						<< "+" << TARNIndianColNamesTable[TARNIndianColAlias::Coupon] << "[i]";
	rowDescVec[TARNIndianColAlias::SumCoupon] = sumCouponDesc.str();
	rowTypeVec[TARNIndianColAlias::SumCoupon] = ARM_STRING;

	CC_Ostringstream	isAliveDesc;
	if (isFirstCpnEvent)
		isAliveDesc << "1";
	else
		isAliveDesc << "if(" << TARNIndianColNamesTable[TARNIndianColAlias::SumCoupon] << "[i]"
					<< "+" << TARNIndianColNamesTable[TARNIndianColAlias::Fees] << "[i]"
					<< "<" << TARNIndianColNamesTable[TARNIndianColAlias::Target] << "[i],1,0)";
	rowDescVec[TARNIndianColAlias::IsAlive] = isAliveDesc.str();
	rowTypeVec[TARNIndianColAlias::IsAlive] = ARM_STRING;

	CC_Ostringstream	paidCouponDesc;
	paidCouponDesc	<< "(" << TARNIndianColNamesTable[TARNIndianColAlias::TwiceCallPut] << "[i]"
					<< "-" << TARNIndianColNamesTable[TARNIndianColAlias::Coupon] << "[i])"
					<< "*DF(" << YC_BASIS_KEY_NAME << itsForeignCcy[0].GetCcyName() << ","
					<< TARNIndianColNamesTable[TARNIndianColAlias::PayDate] << "[i])"
					<< "*" << TARNIndianColNamesTable[TARNIndianColAlias::Notional] << "[i]";
	rowDescVec[TARNIndianColAlias::PaidCoupon] = paidCouponDesc.str();
	rowTypeVec[TARNIndianColAlias::PaidCoupon] = ARM_STRING;

	CC_Ostringstream	realCouponDesc;
	if (isFirstCpnEvent)
		realCouponDesc	<< TARNIndianColNamesTable[TARNIndianColAlias::PaidCoupon] << "[i]";
	else
		realCouponDesc	<< TARNIndianColNamesTable[TARNIndianColAlias::IsAlive] << "[i]"
						<< "*" << TARNIndianColNamesTable[TARNIndianColAlias::PaidCoupon] << "[i]";
	rowDescVec[TARNIndianColAlias::RealCoupon] = realCouponDesc.str();
	rowTypeVec[TARNIndianColAlias::RealCoupon] = ARM_STRING;

	CC_Ostringstream	probaDesc;
	probaDesc	<< "UNPAY(1-" << TARNIndianColNamesTable[TARNIndianColAlias::IsAlive] << "[i]" << ")";
	rowDescVec[TARNIndianColAlias::Proba] = probaDesc.str();
	rowTypeVec[TARNIndianColAlias::Proba] = ARM_STRING;

	CC_Ostringstream	tarnDesc;
	tarnDesc	<< TARNIndianColNamesTable[TARNIndianColAlias::RealCoupon] << "[i]";
	rowDescVec[TARNIndianColAlias::Proba] = tarnDesc.str();
	rowTypeVec[TARNIndianColAlias::Proba] = ARM_STRING;

	return ARM_RowInfo(rowDescVec, rowTypeVec);
}

////////////////////////////////////////////////////
///	Class   : ARM_TARNCalculatorIndian
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_TARNCalculatorIndian::toString(const string& indent, const string& nextIndent) const
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

