/*!
 *
 * Copyright (c) IXIS CIB August 2006 Paris
 *
 *	\file tarnfxcalculator.cpp
 *
 *  \brief file for the PRDKO Switcher
 *	\author  P. Lam
 *	\version 1.0
 *	\date June 2007
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpcalculators/prdkoswitchercalculator.h"
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
////#include "gpbase/datestripconvert.h"
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


const string ARM_PRDKOSwitcherCalculator::SwitcherColNamesTable[] =
{
	"ResetDate",
	"PayDate",
	"IT",
	"Barrier",
	"CallBarrier",
	"Leverage",
	"Notional",
	"SpotFwd",
	"IsAlive",
	"Coupon",
	"RealCoupon",
	"Proba",
	"Option"
};


const int ARM_PRDKOSwitcherCalculator::ProductToPriceColumns[] =
{
	ARM_PRDKOSwitcherCalculator::SpotFwd,
	ARM_PRDKOSwitcherCalculator::IsAlive,
	ARM_PRDKOSwitcherCalculator::Coupon,
	ARM_PRDKOSwitcherCalculator::RealCoupon,
	ARM_PRDKOSwitcherCalculator::Proba,
	ARM_PRDKOSwitcherCalculator::Option
};


/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDKOSwitcherCalculator
///	Routine: Copy constructor
///	Returns: a built object
///	Action : builds the object
/////////////////////////////////////////////////////////////////
ARM_PRDKOSwitcherCalculator::ARM_PRDKOSwitcherCalculator( const ARM_PRDKOSwitcherCalculator& rhs )
:ARM_HybridIRFXCalculator(rhs)
{
	itsBarrier		= rhs.itsBarrier;
	itsCallBarrier	= rhs.itsCallBarrier;
	itsLeverage		= rhs.itsLeverage;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDKOSwitcherCalculator
///	Routine: Destructor
///	Returns: void
///	Action : destroys the object
/////////////////////////////////////////////////////////////////
ARM_PRDKOSwitcherCalculator::~ARM_PRDKOSwitcherCalculator()
{
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDKOSwitcherCalculator
///	Routine: Constructor
///	Returns: void
///	Action : builds the object (stand alone version)
/////////////////////////////////////////////////////////////////
ARM_PRDKOSwitcherCalculator::ARM_PRDKOSwitcherCalculator(const ARM_Date& asOfDate,
														 const ARM_Date& startDate,
														 const ARM_Date& endDate,
														 const ARM_Currency& DomCcy,
														 const vector<ARM_Currency>& ForCcy,
														 const ARM_Currency& CpnCcy,
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
														 const ARM_Curve& barrier,
														 const ARM_Curve& callBarrier,
														 const ARM_Curve& leverage,
														 int intermediatePrices,
														 const ARM_StringVector& columnsToPrice,
														 ARM_FixingSched* pastFixings)
:
	ARM_HybridIRFXCalculator(asOfDate,
							 startDate,
							 endDate,
							 DomCcy,
							 ForCcy,
							 CpnCcy,
							 NULL,
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
							 NULL,
							 NULL,
							 "",
							 intermediatePrices,
							 columnsToPrice,
							 ARM_TARNFXPayoffType::SWITCHER,
							 pastFixings),
	itsBarrier(barrier),
	itsCallBarrier(callBarrier),
	itsLeverage(leverage)
{
	/// Create the Generic Security
	CreateAndSetDealDescription(itsPayModelName, itsColumnsToPrice, ARM_CstManagerPtr(), false);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDKOSwitcherCalculator
///	Routine: CheckData
///	Returns: void
///	Action : check if TARN FX data are consistent
/////////////////////////////////////////////////////////////////
void ARM_PRDKOSwitcherCalculator::CheckData()
{
	if (itsStartDate >= itsEndDate)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : The start date frequency should be before the end date.");
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDKOSwitcherCalculator
///	Routine: CheckMktData
///	Returns: void
///	Action : check if TARN FX market data are consistent
/////////////////////////////////////////////////////////////////
void ARM_PRDKOSwitcherCalculator::CheckMktData()
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

//	if (IsNormal())
		ForexKey = FOREX_KEY_NAME + itsForeignCcy[i].GetCcyName() + "/" + itsDomesticCcy.GetCcyName();
//	else
//		ForexKey = FOREX_KEY_NAME + itsDomesticCcy.GetCcyName() + "/" + itsForeignCcy[i].GetCcyName();
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

//	if (IsNormal())
		FxModelKey = FXMODEL_KEY_NAME + itsForeignCcy[i].GetCcyName() + "/" + itsDomesticCcy.GetCcyName();
//	else
//		FxModelKey = FXMODEL_KEY_NAME + itsDomesticCcy.GetCcyName() + "/" + itsForeignCcy[i].GetCcyName();
	ARM_MixtureModel_Fx* fxMixModel = dynamic_cast< ARM_MixtureModel_Fx* >( GetMktDataManager()->GetData(FxModelKey) );
	ARM_BSModel* fxBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(FxModelKey) );

	if (!fxMixModel && (itsModelType==Model1IRFX))
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : 1IRFX model needs fx Mixture model for key.");

	if(!fxMixModel && !fxBSModel)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : fx Mixture or BS model for key=" + FxModelKey + " is expected in the Market Data Manager");
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDKOSwitcherCalculator
///	Routine: ColumnNames
///	Returns: ARM_RowInfo
///	Action : create the column names of the deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_PRDKOSwitcherCalculator::ColumnNames() const
{
    size_t	colNamesSize = sizeof(SwitcherColNamesTable)/sizeof(SwitcherColNamesTable[0]);
    vector<string>				colNamesVec(colNamesSize);
    vector<ARM_GP_VALUE_TYPE>	colTypeVec(colNamesSize, ARM_STRING); 

    for(size_t i=0; i<colNamesSize; ++i)
        colNamesVec[i] = SwitcherColNamesTable[i];

    ARM_RowInfo	rowInfo(colNamesVec, colTypeVec);

    return rowInfo;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDKOSwitcherCalculator
///	Routine: InitPriceableColumns
///	Returns: void
///	Action : initialise to 0 column to be able to sum
/////////////////////////////////////////////////////////////////
void ARM_PRDKOSwitcherCalculator::InitPriceableColumns(vector<string>& rowDescVec, vector<ARM_GP_VALUE_TYPE>& rowTypeVec) const
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
///	Class  : ARM_PRDKOSwitcherCalculator
///	Routine: DatesStructure
///	Returns: ARM_DateStripCombiner
///	Action : create the list of all event/notice dates of the TARN FX.
//////////////////////////////////////////////////////////////////////
ARM_DateStripCombiner ARM_PRDKOSwitcherCalculator::DatesStructure() const
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

	//skip last flow !!!
	itsCpnDateStrip->ResizeAndBuilt(0, itsCpnDateStrip->size()-1); 

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

	const_cast<ARM_PRDKOSwitcherCalculator*>(this)->itsFirstEventIdx = firstEventIdx;

    return	EventSchedule;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDKOSwitcherCalculator
///	Routine: MiddleRows
///	Returns: ARM_RowInfo
///	Action : create a row of a deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_PRDKOSwitcherCalculator::MiddleRows( size_t eventIdx, const ARM_DateStripCombiner& datesStructure ) const
{
	ARM_DateStripPtr	vDateStrip = datesStructure.GetDateStrip(0);
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	size_t	eventSize = vDateStrip->GetResetDates()->size();
	size_t	descSize = sizeof(SwitcherColNamesTable)/sizeof(SwitcherColNamesTable[0]);

	vector<string>				rowDescVec(descSize);
	vector<ARM_GP_VALUE_TYPE>	rowTypeVec(descSize, ARM_MISSING_TYPE); 

    /// Set default 0 value for each column to be able to sum it
    InitPriceableColumns(rowDescVec, rowTypeVec);
	int i = 0;

	bool isFirstEvent = ((itsCpnTiming == K_ARREARS)&&(itsFirstEventIdx==itsNbPastFixings)?(eventIdx == itsNbPastFixings+1):(eventIdx == itsFirstEventIdx));

    CC_Ostringstream	resetDateDesc;
	double vResetDate = (*(vDateStrip->GetResetDates()))[eventIdx];
    resetDateDesc << CC_NS(std, fixed) << vResetDate;
    rowDescVec[SwitcherColAlias::ResetDate] = resetDateDesc.str();
    rowTypeVec[SwitcherColAlias::ResetDate] = ARM_DATE_TYPE;

	CC_Ostringstream	payDateDesc;
	double vPayDate = (*(vDateStrip->GetFlowStartDates()))[eventIdx]; // WARNING !
	payDateDesc << CC_NS(std, fixed) << vPayDate;
	rowDescVec[SwitcherColAlias::PayDate] = payDateDesc.str();
	rowTypeVec[SwitcherColAlias::PayDate] = ARM_DATE_TYPE;

	CC_Ostringstream	ITDesc;
	ITDesc << CC_NS(std, fixed) << (*(vDateStrip->GetInterestTerms()))[eventIdx];
	rowDescVec[SwitcherColAlias::IT] = ITDesc.str();
	rowTypeVec[SwitcherColAlias::IT] = ARM_DOUBLE;

	CC_Ostringstream	barrierDesc;
	barrierDesc << CC_NS(std, fixed) << itsBarrier.Interpolate(vResetDate-asOfDate);
	rowDescVec[SwitcherColAlias::Barrier] = barrierDesc.str();
	rowTypeVec[SwitcherColAlias::Barrier] = ARM_DOUBLE;

	CC_Ostringstream	callBarrierDesc;
	callBarrierDesc << CC_NS(std, fixed) << itsCallBarrier.Interpolate(vResetDate-asOfDate);
	rowDescVec[SwitcherColAlias::CallBarrier] = callBarrierDesc.str();
	rowTypeVec[SwitcherColAlias::CallBarrier] = ARM_DOUBLE;
	
	CC_Ostringstream	leverageDesc;
	leverageDesc << CC_NS(std, fixed) << itsLeverage.Interpolate(vResetDate-asOfDate);
	rowDescVec[SwitcherColAlias::Leverage] = leverageDesc.str();
	rowTypeVec[SwitcherColAlias::Leverage] = ARM_DOUBLE;

	CC_Ostringstream	notionalDesc;
	notionalDesc << CC_NS(std, fixed) << itsCpnNominalCv.Interpolate(vPayDate-asOfDate);
	rowDescVec[SwitcherColAlias::Notional] = notionalDesc.str();
	rowTypeVec[SwitcherColAlias::Notional] = ARM_DOUBLE;

	CC_Ostringstream	spotFwdDesc;
	spotFwdDesc << CC_NS(std, fixed) << const_cast<ARM_PRDKOSwitcherCalculator*>(this)->FXValue(vDateStrip, eventIdx);
	rowDescVec[SwitcherColAlias::SpotFwd] = spotFwdDesc.str();
	rowTypeVec[SwitcherColAlias::SpotFwd] = ARM_STRING;

	CC_Ostringstream	isAliveDesc;
	isAliveDesc << "if(" << SwitcherColNamesTable[SwitcherColAlias::SpotFwd] << "[i]"
				<< "<" << SwitcherColNamesTable[SwitcherColAlias::CallBarrier] << "[i],1,0)";
	if (!isFirstEvent)
		isAliveDesc << "*" << SwitcherColNamesTable[SwitcherColAlias::IsAlive] << "[i-1]";

	rowDescVec[SwitcherColAlias::IsAlive] = isAliveDesc.str();
	rowTypeVec[SwitcherColAlias::IsAlive] = ARM_STRING;

	CC_Ostringstream	couponDesc;
	couponDesc	<< SwitcherColNamesTable[SwitcherColAlias::Leverage] << "[i]"
				<< "*Max(" << SwitcherColNamesTable[SwitcherColAlias::Barrier] << "[i]"
				<< "-" << SwitcherColNamesTable[SwitcherColAlias::SpotFwd] << "[i],0.0)"
				<< "*" << SwitcherColNamesTable[SwitcherColAlias::IT] << "[i]"
				<< "*DF(" << YC_KEY_NAME << itsCouponCcy.GetCcyName() << ","
				<< SwitcherColNamesTable[SwitcherColAlias::PayDate] << "[i])"
				<< "*" << SwitcherColNamesTable[SwitcherColAlias::Notional] << "[i]";
	rowDescVec[SwitcherColAlias::Coupon] = couponDesc.str();
	rowTypeVec[SwitcherColAlias::Coupon] = ARM_STRING;

	CC_Ostringstream	realCouponDesc;
	realCouponDesc	<< SwitcherColNamesTable[SwitcherColAlias::IsAlive] << "[i]"
					<< "*" << SwitcherColNamesTable[SwitcherColAlias::Coupon] << "[i]";
	rowDescVec[SwitcherColAlias::RealCoupon] = realCouponDesc.str();
	rowTypeVec[SwitcherColAlias::RealCoupon] = ARM_STRING;

	CC_Ostringstream	probaDesc;
	probaDesc	<< "UNPAY(" << SwitcherColNamesTable[SwitcherColAlias::IsAlive] << "[i]" << ")";
	rowDescVec[SwitcherColAlias::Proba] = probaDesc.str();
	rowTypeVec[SwitcherColAlias::Proba] = ARM_STRING;

	CC_Ostringstream	optionDesc;
	optionDesc	<< SwitcherColNamesTable[SwitcherColAlias::RealCoupon] << "[i]"
				<< "-" << SwitcherColNamesTable[SwitcherColAlias::Coupon] << "[i]";
	rowDescVec[SwitcherColAlias::Option] = optionDesc.str();
	rowTypeVec[SwitcherColAlias::Option] = ARM_STRING;

	return ARM_RowInfo(rowDescVec, rowTypeVec);
}

////////////////////////////////////////////////////
///	Class   : ARM_PRDKOSwitcherCalculator
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_PRDKOSwitcherCalculator::toString(const string& indent, const string& nextIndent) const
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