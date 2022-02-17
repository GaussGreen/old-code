/*!
 *
 * Copyright (c) IXIS CIB August 2006 Paris
 *
 *	\file tarnfxcalculator.cpp
 *
 *  \brief file for the TARN Calculator Indian
 *	\author  P. Lam
 *	\version 1.0
 *	\date May 2007
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


/// Default MDM key names
const string YC_KEY_NAME		= "YC_";
const string YC_BASIS_KEY_NAME	= "YC_BASIS_";
const string FOREX_KEY_NAME		= "FOREX_";
const string OSWMODEL_KEY_NAME	= "OSWMOD_";
const string FXMODEL_KEY_NAME	= "FXMOD_";
const string CORREL_KEY_NAME	= "CORREL_";
const string MRS_KEY_NAME		= "MRS_";
const string QFX_KEY_NAME		= "Q_";


const string ARM_TARNCalculatorIndian::IndianColNamesTable[] =
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
	"Short",
	"Long",
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
	ARM_TARNCalculatorIndian::Short,
	ARM_TARNCalculatorIndian::Long,
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
:ARM_HybridIRFXCalculator(rhs)
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
ARM_TARNCalculatorIndian::ARM_TARNCalculatorIndian (const ARM_Date& asOfDate,
													const ARM_Date& startDate,
													const ARM_Date& endDate,
													const ARM_Currency& CpnCcy,
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
													const ARM_Curve& strike,
													const ARM_Curve& target,
													double epsilon,
													const ARM_Curve& fees,
													int intermediatePrices,
													const ARM_StringVector& columnsToPrice,
													const int optionType,
													const int indianType,
													ARM_FixingSched* pastFixings)
:
	ARM_HybridIRFXCalculator(asOfDate,
							 startDate,
							 endDate,
							 CpnCcy,
							 ForCcy,
							 &CpnCcy,
							 payRec,
							 cpnDayCount,
							 cpnFreq,
							 cpnResetGap,
							 cpnResetCal,
							 cpnPayCal,
							 stubRule,
							 K_ADVANCE, //reset timing not relevant : always ADV
							 cpnIntRule,
							 &cpnNominal,
							 NULL,
							 NULL,
							 NULL,
							 NULL,
							 NULL,
							 cpnFreq,
							 KNOBASE,
							 NULL,
							 NULL,
							 0,
							 0,
							 0.0,
							 target,
							 &fees,
							 "",
							 intermediatePrices,
							 columnsToPrice,
							 ARM_TARNFXPayoffType::INDIANFX,
							 pastFixings),
	itsBarrierUp(barrierUp),
	itsEpsilon(epsilon),
	itsStrike(strike),
	itsOptionType(optionType),
	itsIndianType(indianType)
{
	// barrier down can be null
	itsBarrierDown = (&barrierDown) ? barrierDown : ARM_Curve();

	/// Create the Generic Security
	CreateAndSetDealDescription(itsPayModelName, itsColumnsToPrice, ARM_CstManagerPtr(), false);
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
    size_t	colNamesSize = sizeof(IndianColNamesTable)/sizeof(IndianColNamesTable[0]);
    vector<string>				colNamesVec(colNamesSize);
    vector<ARM_GP_VALUE_TYPE>	colTypeVec(colNamesSize, ARM_STRING); 

    for(size_t i=0; i<colNamesSize; ++i)
        colNamesVec[i] = IndianColNamesTable[i];

    ARM_RowInfo	rowInfo(colNamesVec, colTypeVec);

    return rowInfo;
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

	int i = 0;
	for (; i<(sizeof(ProductToPriceColumns)/sizeof(ProductToPriceColumns[0])); i++)
	{
		rowDescVec[ProductToPriceColumns[i]] = zeroValue;
		rowTypeVec[ProductToPriceColumns[i]] = ARM_DOUBLE;
	}
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
    SchedVect[CPNFX_SCHED] = &*itsCpnDateStrip;

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
	size_t	descSize = sizeof(IndianColNamesTable)/sizeof(IndianColNamesTable[0]);

	vector<string>				rowDescVec(descSize);
	vector<ARM_GP_VALUE_TYPE>	rowTypeVec(descSize, ARM_MISSING_TYPE); 

    /// Set default 0 value for each column to be able to sum it
    InitPriceableColumns(rowDescVec, rowTypeVec);
	int i = 0;

	bool isFirstEvent = ((itsCpnTiming == K_ARREARS)&&(itsFirstEventIdx==itsNbPastFixings)?(eventIdx == itsNbPastFixings+1):(eventIdx == itsFirstEventIdx));

    CC_Ostringstream	resetDateDesc;
	double vResetDate = (*(vDateStrip->GetResetDates()))[eventIdx];
    resetDateDesc << CC_NS(std, fixed) << vResetDate;
    rowDescVec[IndianColAlias::ResetDate] = resetDateDesc.str();
    rowTypeVec[IndianColAlias::ResetDate] = ARM_DATE_TYPE;

	CC_Ostringstream	payDateDesc;
	double vPayDate = (*(vDateStrip->GetFlowStartDates()))[eventIdx]; // WARNING !
	payDateDesc << CC_NS(std, fixed) << vPayDate;
	rowDescVec[IndianColAlias::PayDate] = payDateDesc.str();
	rowTypeVec[IndianColAlias::PayDate] = ARM_DATE_TYPE;

	CC_Ostringstream	barrierUpDesc;
	barrierUpDesc << CC_NS(std, fixed) << itsBarrierUp.Interpolate(vResetDate-asOfDate);
	rowDescVec[IndianColAlias::BarrierUp] = barrierUpDesc.str();
	rowTypeVec[IndianColAlias::BarrierUp] = ARM_DOUBLE;


	if (itsIndianType == ARM_TARNCalculatorIndian::DownUp)
	{
		CC_Ostringstream	barrierDownDesc;
		barrierDownDesc << CC_NS(std, fixed) << itsBarrierDown.Interpolate(vResetDate-asOfDate);
		rowDescVec[IndianColAlias::BarrierDown] = barrierDownDesc.str();
		rowTypeVec[IndianColAlias::BarrierDown] = ARM_DOUBLE;
	}

	CC_Ostringstream	strikeDesc;
	strikeDesc << CC_NS(std, fixed) << itsStrike.Interpolate(vResetDate-asOfDate);
	rowDescVec[IndianColAlias::Strike] = strikeDesc.str();
	rowTypeVec[IndianColAlias::Strike] = ARM_DOUBLE;

	CC_Ostringstream	targetRedemptionDesc; 
	targetRedemptionDesc << CC_NS(std, fixed) << itsTarget.Interpolate(vResetDate-asOfDate);
	rowDescVec[IndianColAlias::Target] = targetRedemptionDesc.str();
	rowTypeVec[IndianColAlias::Target] = ARM_DOUBLE;

	CC_Ostringstream	notionalDesc;
	notionalDesc << CC_NS(std, fixed) << itsCpnNominalCv.Interpolate(vPayDate-asOfDate);
	rowDescVec[IndianColAlias::Notional] = notionalDesc.str();
	rowTypeVec[IndianColAlias::Notional] = ARM_DOUBLE;

	CC_Ostringstream	feesDesc;
	feesDesc << CC_NS(std, fixed) << itsFees.Interpolate(vPayDate-asOfDate);
	rowDescVec[IndianColAlias::Fees] = feesDesc.str();
	rowTypeVec[IndianColAlias::Fees] = ARM_DOUBLE;

	CC_Ostringstream	spotFwdDesc;
	spotFwdDesc << CC_NS(std, fixed) << const_cast<ARM_TARNCalculatorIndian*>(this)->FXValue(vDateStrip, eventIdx);
	rowDescVec[IndianColAlias::SpotFwd] = spotFwdDesc.str();
	rowTypeVec[IndianColAlias::SpotFwd] = ARM_STRING;

	if (itsIndianType == DownUp)
	{
		CC_Ostringstream	twiceCallPutDesc;
		twiceCallPutDesc	<< "2*(Max(" << IndianColNamesTable[IndianColAlias::SpotFwd] << "[i]"
							<< "-" << IndianColNamesTable[IndianColAlias::BarrierUp] << "[i],0.0)"
							<< "+Max(" << IndianColNamesTable[IndianColAlias::BarrierDown] << "[i]"
							<< "-" << IndianColNamesTable[IndianColAlias::SpotFwd] << "[i],0.0))";
		rowDescVec[IndianColAlias::TwiceCallPut] = twiceCallPutDesc.str();
		rowTypeVec[IndianColAlias::TwiceCallPut] = ARM_STRING;

		CC_Ostringstream	callMinusCallDesc;
		callMinusCallDesc	<< "Max(" << IndianColNamesTable[IndianColAlias::SpotFwd] << "[i]"
							<< "-" << IndianColNamesTable[IndianColAlias::BarrierDown] << "[i],0.0)"
							<< "-Max(" << IndianColNamesTable[IndianColAlias::SpotFwd] << "[i]"
							<< "-" << IndianColNamesTable[IndianColAlias::Strike] << "[i],0.0)";
		rowDescVec[IndianColAlias::CallMinusCall] = callMinusCallDesc.str();
		rowTypeVec[IndianColAlias::CallMinusCall] = ARM_STRING;

		CC_Ostringstream	putMinusPutDesc;
		putMinusPutDesc	<< "Max(" << IndianColNamesTable[IndianColAlias::BarrierUp] << "[i]"
						<< "-" << IndianColNamesTable[IndianColAlias::SpotFwd] << "[i],0.0)"
						<< "-Max(" << IndianColNamesTable[IndianColAlias::Strike] << "[i]"
						<< "-" << IndianColNamesTable[IndianColAlias::SpotFwd] << "[i],0.0)";
		rowDescVec[IndianColAlias::PutMinusPut] = putMinusPutDesc.str();
		rowTypeVec[IndianColAlias::PutMinusPut] = ARM_STRING;

		CC_Ostringstream	epsilonDesc;
		epsilonDesc	<< "0.01*" << IndianColNamesTable[IndianColAlias::Strike] << "[i]";
		rowDescVec[IndianColAlias::Epsilon] = epsilonDesc.str();
		rowTypeVec[IndianColAlias::Epsilon] = ARM_STRING;

		CC_Ostringstream	callDigitalDesc;
		callDigitalDesc	<< "(" << IndianColNamesTable[IndianColAlias::Strike] << "[i]"
						<< "-" << IndianColNamesTable[IndianColAlias::BarrierDown] << "[i])"
						<< "*(Max(" << IndianColNamesTable[IndianColAlias::SpotFwd] << "[i]"
						<< "-(" << IndianColNamesTable[IndianColAlias::Strike] << "[i]"
						<< "-" << IndianColNamesTable[IndianColAlias::Epsilon] << "[i]),0.0)"
						<< "-Max(" << IndianColNamesTable[IndianColAlias::SpotFwd] << "[i]"
						<< "-(" << IndianColNamesTable[IndianColAlias::Strike] << "[i]"
						<< "+" << IndianColNamesTable[IndianColAlias::Epsilon] << "[i]),0.0))"
						<< "/(2*" << IndianColNamesTable[IndianColAlias::Epsilon] << "[i])";
		rowDescVec[IndianColAlias::CallDigital] = callDigitalDesc.str();
		rowTypeVec[IndianColAlias::CallDigital] = ARM_STRING;

		CC_Ostringstream	putDigitalDesc;
		putDigitalDesc	<< "(" << IndianColNamesTable[IndianColAlias::BarrierUp] << "[i]"
						<< "-" << IndianColNamesTable[IndianColAlias::Strike] << "[i])"
						<< "*(Max((" << IndianColNamesTable[IndianColAlias::Strike] << "[i]"
						<< "+" << IndianColNamesTable[IndianColAlias::Epsilon] << "[i])"
						<< "-" << IndianColNamesTable[IndianColAlias::SpotFwd] << "[i],0.0)"
						<< "-Max((" << IndianColNamesTable[IndianColAlias::Strike] << "[i]"
						<< "-" << IndianColNamesTable[IndianColAlias::Epsilon] << "[i])"
						<< "-" << IndianColNamesTable[IndianColAlias::SpotFwd] << "[i],0.0))"
						<< "/(2*" << IndianColNamesTable[IndianColAlias::Epsilon] << "[i])";
		rowDescVec[IndianColAlias::PutDigital] = putDigitalDesc.str();
		rowTypeVec[IndianColAlias::PutDigital] = ARM_STRING;

		CC_Ostringstream	couponDesc;
		couponDesc	<< IndianColNamesTable[IndianColAlias::CallMinusCall] << "[i]"
					<< "-" << IndianColNamesTable[IndianColAlias::CallDigital] << "[i]"
					<< "+" << IndianColNamesTable[IndianColAlias::PutMinusPut] << "[i]"
					<< "-" << IndianColNamesTable[IndianColAlias::PutDigital] << "[i]";
		rowDescVec[IndianColAlias::Coupon] = couponDesc.str();
		rowTypeVec[IndianColAlias::Coupon] = ARM_STRING;

		CC_Ostringstream	sumCouponDesc;
		sumCouponDesc << IndianColNamesTable[IndianColAlias::Coupon] << "[i]";
		if (!isFirstEvent)
			sumCouponDesc	<< "+" << IndianColNamesTable[IndianColAlias::SumCoupon] << "[i-1]";

		rowDescVec[IndianColAlias::SumCoupon] = sumCouponDesc.str();
		rowTypeVec[IndianColAlias::SumCoupon] = ARM_STRING;

		CC_Ostringstream	paidCouponDesc;
		paidCouponDesc	<< "(" << IndianColNamesTable[IndianColAlias::TwiceCallPut] << "[i]"
						<< "-" << IndianColNamesTable[IndianColAlias::Coupon] << "[i])"
						<< "*DF(" << YC_BASIS_KEY_NAME << itsDomesticCcy.GetCcyName() << ","
						<< IndianColNamesTable[IndianColAlias::PayDate] << "[i])"
						<< "*" << IndianColNamesTable[IndianColAlias::Notional] << "[i]";
		rowDescVec[IndianColAlias::PaidCoupon] = paidCouponDesc.str();
		rowTypeVec[IndianColAlias::PaidCoupon] = ARM_STRING;
	}
	else
	{
		CC_Ostringstream	shortDesc;
		if (itsOptionType == K_PUT)
			shortDesc	<< "Max(" << IndianColNamesTable[IndianColAlias::Strike] << "[i]"
						<< "-" << IndianColNamesTable[IndianColAlias::SpotFwd] << "[i],0.0)";
		else
			shortDesc	<< "Max(" << IndianColNamesTable[IndianColAlias::SpotFwd] << "[i]"
						<< "-" << IndianColNamesTable[IndianColAlias::Strike] << "[i],0.0)";
		rowDescVec[IndianColAlias::Short] = shortDesc.str();
		rowTypeVec[IndianColAlias::Short] = ARM_STRING;

		CC_Ostringstream	longDesc;
		longDesc	<< "2*if(" << IndianColNamesTable[IndianColAlias::SpotFwd] << "[i]"
					<< ">" << IndianColNamesTable[IndianColAlias::BarrierUp] << "[i],"
					<< IndianColNamesTable[IndianColAlias::SpotFwd] << "[i]"
					<< "-" << IndianColNamesTable[IndianColAlias::Strike] << "[i],0.0)";
		rowDescVec[IndianColAlias::Long] = longDesc.str();
		rowTypeVec[IndianColAlias::Long] = ARM_STRING;

		CC_Ostringstream	sumCouponDesc;
		if (itsIndianType == Trigger)
		{
			sumCouponDesc << IndianColNamesTable[IndianColAlias::Short] << "[i]";
		}
		else
		{
			sumCouponDesc	<< "if(" << IndianColNamesTable[IndianColAlias::SpotFwd] << "[i]"
							<< "<=Strike[i],1,0)";
		}

		if (!isFirstEvent)
			sumCouponDesc	<< "+" << IndianColNamesTable[IndianColAlias::SumCoupon] << "[i-1]";
		rowDescVec[IndianColAlias::SumCoupon] = sumCouponDesc.str();
		rowTypeVec[IndianColAlias::SumCoupon] = ARM_STRING;

		CC_Ostringstream	paidCouponDesc;
		paidCouponDesc	<< "(Long[i]-Short[i])"
						<< "*DF(" << YC_BASIS_KEY_NAME << itsDomesticCcy.GetCcyName() << ","
						<< IndianColNamesTable[IndianColAlias::PayDate] << "[i])"
						<< "*" << IndianColNamesTable[IndianColAlias::Notional] << "[i]";
		rowDescVec[IndianColAlias::PaidCoupon] = paidCouponDesc.str();
		rowTypeVec[IndianColAlias::PaidCoupon] = ARM_STRING;
	}

	CC_Ostringstream	isAliveDesc;
	if (itsIndianType == DownUp)
	{
		if (isFirstEvent)
			isAliveDesc << "1";
		else
			isAliveDesc << "if(" << IndianColNamesTable[IndianColAlias::SumCoupon] << "[i]"
						<< "+" << IndianColNamesTable[IndianColAlias::Fees] << "[i]"
						<< "<" << IndianColNamesTable[IndianColAlias::Target] << "[i],1,0)";
	}
	else
	{
		isAliveDesc << "if(" << IndianColNamesTable[IndianColAlias::SumCoupon] << "[i]"
					<< "+" << IndianColNamesTable[IndianColAlias::Fees] << "[i]"
					<< "<" << IndianColNamesTable[IndianColAlias::Target] << "[i],1,0)";

		if (!isFirstEvent)
			isAliveDesc << "*" << IndianColNamesTable[IndianColAlias::IsAlive] << "[i-1]";
	}
	rowDescVec[IndianColAlias::IsAlive] = isAliveDesc.str();
	rowTypeVec[IndianColAlias::IsAlive] = ARM_STRING;

	CC_Ostringstream	realCouponDesc;
	realCouponDesc	<< IndianColNamesTable[IndianColAlias::PaidCoupon] << "[i]";
	if (!isFirstEvent)
	{
		if (itsIndianType == DownUp)
			realCouponDesc	<< "*" << IndianColNamesTable[IndianColAlias::IsAlive] << "[i]";
		else
			realCouponDesc	<< "*" << IndianColNamesTable[IndianColAlias::IsAlive] << "[i-1]";
	}
	rowDescVec[IndianColAlias::RealCoupon] = realCouponDesc.str();
	rowTypeVec[IndianColAlias::RealCoupon] = ARM_STRING;

	CC_Ostringstream	probaDesc;
	probaDesc	<< "UNPAY(1-" << IndianColNamesTable[IndianColAlias::IsAlive] << "[i]" << ")";
	rowDescVec[IndianColAlias::Proba] = probaDesc.str();
	rowTypeVec[IndianColAlias::Proba] = ARM_STRING;

	CC_Ostringstream	tarnDesc;
	tarnDesc	<< IndianColNamesTable[IndianColAlias::RealCoupon] << "[i]";
	rowDescVec[IndianColAlias::TARN] = tarnDesc.str();
	rowTypeVec[IndianColAlias::TARN] = ARM_STRING;

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