/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ccscalculator.cpp
 *
 *  \brief file for the Callable cross currency swap Calculator
 *	\author  K.Belkheir & E.Ezzine
 *	\version 1.0
 *	\date January 2007
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpcalculators/ccscalculator.h"

/// gpbase
#include "gpbase/env.h"
#include "gpbase/autocleaner.h"
#include "gpbase/ostringstream.h"
#include "gpbase/datestripcombiner.h"
#include "gpbase/datestrip.h"
#include "gpbase/singleton.h"
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
#include "gpinfra/pricerinfo.h"
#include "gpinfra/cstmanager.h"
#include "gpinfra/discretisationscheme.h"

/// gpcalib
#include "gpcalib/calibmethod.h"
#include "gpcalib/modelfitter.h"

/// gpmodels
#include "gpmodels/2irfxModel.h"
#include "gpmodels/q1f.h"
#include "gpmodels/q1f_fx.h"
#include "gpmodels/modelparamsq1f.h"
#include "gpmodels/eqfx_modelfactory.h"
#include "gpmodels/ModelParams_EqFxBase.h"
#include "gpmodels/forwardmarginbasis.h"

/// gpnummethods
#include "gpnummethods/treebase.h"
#include "gpnummethods/treefactory.h"
#include "gpnummethods/scheduler.h"
#include "gpnummethods/sampler.h"
#include "gpnummethods/truncator.h"
#include "gpnummethods/reconnector.h"
#include "gpnummethods/smoother.h"


/// kernel
#include <inst/forex.h>
#include <inst/swaption.h>
#include <inst/option.h>
#include <inst/portfolio.h>

/// STL
#include <iomanip> /// for setprecision()
#include <list>
CC_USING_NS(std,list)


CC_BEGIN_NAMESPACE( ARM )

const double NON_CALL_FEE			= 1.0e15;

/// H&W vol range [10bp,500bp]
const double HWVOL_LOWER_BOUND      = 0.001;
const double HWVOL_UPPER_BOUND      = 0.05;

/// Q vol range [2.5%,100%]
const double QVOL_LOWER_BOUND       = 0.025;
const double QVOL_UPPER_BOUND       = 1.0;

/// 10-3 bp of vega to be selected in portfolio for IR volatility bootstrapping 
const double IR_VEGA_MIN_TO_SELECT  = 1.0e-7;
const double OSW_DEFAULT_WEIGHT     = 1.0;
const double OSW_DEFAULT_PRICE      = 1.0e+10;

/// 10-3 bp of vega to be selected in portfolio for FX volatility bootstrapping 
const double FX_VEGA_MIN_TO_SELECT  = 1.0e-7;
const double FX_DEFAULT_WEIGHT      = 1.0;
const double FX_DEFAULT_PRICE       = 1.0e+30;

/// For ATS minimum volatility search
const double FX_MIN_MONEYNESS       = 0.5;
const double FX_MAX_MONEYNESS       = 1.75;
const double FX_STEP_MONEYNESS      = 0.05;

/// Default MDM key names
const string YC_KEY_NAME                    = "YC_";
const string YC_BASIS_KEY_NAME              = "YC_BASIS_";
const string FOREX_KEY_NAME                 = "FOREX_";
const string OSWMODEL_KEY_NAME              = "OSWMOD_";
const string FXMODEL_KEY_NAME               = "FXMOD_";
const string CORREL_KEY_NAME                = "CORREL_";
const string MRS_KEY_NAME                   = "MRS_";
const string Q_KEY_NAME                     = "Q_";

const string ARM_CCSCalculator::CCSColNamesTable [] =
{
	"EventDate",		
	"EndDate",
	"Spot",		
	"ForStartDate",
	"ForEndDate",
	"ForeignLeg",
	"ForNominal",	
	"ForExchangeNotional",
	"DomStartDate",
	"DomEndDate",
	"DomesticLeg",
	"DomNominal",	
	"DomExchangeNotional",
	"ForeignFlow",
	"DomesticFlow",
	"CCSSwapFlow",
	"BasisCCSSwap",
	"CCSwap",
	"FirstCCSSwap",        
	"CCSBermuda",
	"CCSOption"
};

const string ARM_CCSCalculator::CCSProfileNamesTable [] =
{
	 "FxDateStripProfile",
	 "DomDateStripProfile",
	 "ForDateStripProfile",
	 "DomNominalProfile",
	 "ForNominalProfile",
     "DomMarginProfile",
	 "ForMarginProfile"
};


/// Maximum calibration errors in % (IR) and real terms (FX)
const double MAX_CALIB_ERR[] = {    0.01,
                                    0.01,
                                    0.0001};

const string CALIB_ERR_MSGE[] = {   ": domestic IR calibration failed",
                                    ": foreign IR calibration failed",
                                    ": forex calibration failed"};


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CCSCalculator
///	Routine: Copy constructor
///	Returns: a built object
///	Action : builds the object
/////////////////////////////////////////////////////////////////
ARM_CCSCalculator::ARM_CCSCalculator( const ARM_CCSCalculator& rhs )
:	ARM_GenCalculator( rhs ),
	
	itsStartDate(rhs.itsStartDate),          
	itsFixedEndDate(rhs.itsFixedEndDate),       
	itsEndDate(rhs.itsEndDate),       
	itsDomesticCcy(rhs.itsDomesticCcy),        
	itsForeignCcy(rhs.itsForeignCcy),		
	itsDomDaycount(rhs.itsDomDaycount),        
	itsDomFreq(rhs.itsDomFreq),            
	itsFxResetGap(rhs.itsFxResetGap),
	itsResetCalendar(rhs.itsResetCalendar),        
	itsPayCalendar(rhs.itsPayCalendar),
	itsForResetCalendar(rhs.itsForResetCalendar),        
	itsForPayCalendar(rhs.itsForPayCalendar),   
	itsStubRule(rhs.itsStubRule),           					
	itsForFreq(rhs.itsForFreq),
	itsForDaycount(rhs.itsForDaycount),
	itsExerciseFreq(rhs.itsExerciseFreq),
	itsNoticeGap(rhs.itsNoticeGap),
	itsPayRec(rhs.itsPayRec),
	itsNbNoCall(rhs.itsNbNoCall),
	itsFees(rhs.itsFees),
	itsDomNominalCv(rhs.itsDomNominalCv),
	itsForNominalCv(rhs.itsForNominalCv),
	itsDomSpreadCv(rhs.itsDomSpreadCv),
	itsForSpreadCv(rhs.itsForSpreadCv),
	itsNbFixFlows(rhs.itsNbFixFlows),   
	itsvDomSpread(rhs.itsvDomSpread),
	itsvForSpread(rhs.itsvForSpread),
	itsvDomNominal(rhs.itsvDomNominal),
	itsvForNominal(rhs.itsvForNominal),
	itsForSize(rhs.itsForSize),
	itsvForIndex(rhs.itsvForIndex),
	itsvFixCpn(rhs.itsvFixCpn),
	itsvCpnIsFixed(rhs.itsvCpnIsFixed),	    
	itsDomSize(rhs.itsDomSize),
	itsvDomIndex(rhs.itsvDomIndex),
	itsvIsExerDate(rhs.itsvIsExerDate),
	itsExerSize(rhs.itsExerSize),
	itsDomDateStrip(rhs.itsDomDateStrip),
	itsForDateStrip(rhs.itsForDateStrip),
	itsForexDateStrip(rhs.itsForexDateStrip),
	itsExerciseDateStrip(rhs.itsExerciseDateStrip),
	itsSchedulerDatas(rhs.itsSchedulerDatas),
	itsTruncatorDatas(rhs.itsTruncatorDatas),
	itsColumnsToPrice(rhs.itsColumnsToPrice),
	itsMarkovianDriftSamplerFlag(rhs.itsMarkovianDriftSamplerFlag),
	itsCalibType(rhs.itsCalibType),
	itsCalibDatas(rhs.itsCalibDatas),
	itsHasBeenComputed(rhs.itsHasBeenComputed)
{
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CCSCalculator
///	Routine: Constructor
///	Returns: void
///	Action : builds the object (stand alone version)
/////////////////////////////////////////////////////////////////
ARM_CCSCalculator::ARM_CCSCalculator(const ARM_Date& asOfDate,
		const ARM_Date& startDate,
		const ARM_Date& fixEndDate,
		const ARM_Date& endDate,
		const ARM_Currency& domCcy,
		const ARM_Currency& forCcy,
		int domDayCount,
		int domFreq,
		const string& domResetCal,
		const string& domPayCal,
		int forFreq,
		int forDayCount,
		const string& forResetCal,
		const string& forPayCal,
		int FxResetGap,
		int stubRule,	
		const ARM_Curve& domMargin,	
		const ARM_Curve& domNominal,		
		const ARM_Curve& forMargin,
		const ARM_Curve& forNominal,				
		int exerciseFreq,
		int noticeGap,
		int payRec,
		size_t nbNoCall,
		const ARM_Curve& fees,
		const ARM_StringVector& productsToPrice)
:	
	ARM_GenCalculator(asOfDate),
	itsStartDate(startDate),          
	itsFixedEndDate(fixEndDate),       
	itsEndDate(endDate),       
	itsDomesticCcy(domCcy),        
	itsForeignCcy(forCcy),		
	itsDomDaycount(domDayCount),        
	itsDomFreq(domFreq),            
	itsFxResetGap(FxResetGap),
	itsResetCalendar(domResetCal),        
	itsPayCalendar(domPayCal),  
	itsForResetCalendar(forResetCal),        
	itsForPayCalendar(forPayCal), 
	itsStubRule(stubRule),           					
	itsForFreq(forFreq),
	itsForDaycount(forDayCount),
	itsExerciseFreq(exerciseFreq),
	itsNoticeGap(noticeGap),
	itsPayRec(payRec),
	itsNbNoCall(nbNoCall),
	itsFees(fees),
	itsDomNominalCv(domNominal),
	itsForNominalCv(forNominal),
	itsDomSpreadCv(domMargin),
	itsForSpreadCv(forMargin),
	itsNbFixFlows(0),   
	itsvDomSpread(0),
	itsvForSpread(0),
	itsvDomNominal(0),
	itsvForNominal(0),
	itsForSize(0),
	itsvForIndex(0),
	itsvFixCpn(0),
	itsvCpnIsFixed(0),	    
	itsDomSize(0),
	itsvDomIndex(0),
	itsvIsExerDate(0),
	itsExerSize(0),
	itsDomDateStrip(0),
	itsForDateStrip(0),
	itsForexDateStrip(0),
	itsExerciseDateStrip(0),
	itsSchedulerDatas(0),
	itsTruncatorDatas(0),
	itsColumnsToPrice(productsToPrice),
	itsMarkovianDriftSamplerFlag(true),
	itsCalibType(ARM_PRCSCalibTypes::ATMCalib),
	itsCalibDatas(0),
	itsHasBeenComputed(false)
{
	DatesStructure();

	/// fix period managing by fixEndDate
	while(itsNbFixFlows < itsDomDateStrip->GetFlowEndDates()->size() && 
		(*itsDomDateStrip->GetFlowEndDates())[itsNbFixFlows] < fixEndDate.GetJulian()+ ARM_GlobalConstant::ARM_SEVENDAYS_LAG )
		itsNbFixFlows++;

	/// Check input datas
    CheckDataAndTimeIt();

    /// Set the domestic=coupon=payment currency
    SetDomesticCcy(const_cast<ARM_Currency&> (domCcy));
	SetForeignCcy(const_cast<ARM_Currency&> (forCcy));

    string domCcyName( domCcy.GetCcyName() );
    string forCcyName( forCcy.GetCcyName() );
    string fxName(forCcyName + "/" + domCcyName);
    string domForFxName("(" + domCcyName + "," + forCcyName + "," + fxName + ")");

    ARM_StringVector mdmKeys(NbKeys);
    mdmKeys[YcDomKey]               = YC_KEY_NAME               + domCcyName;
    mdmKeys[YcForKey]               = YC_KEY_NAME               + forCcyName;
    mdmKeys[ForexKey]               = FOREX_KEY_NAME            + fxName;
    mdmKeys[YcBasisDomKey]          = YC_BASIS_KEY_NAME         + domCcyName;
    mdmKeys[YcBasisForKey]          = YC_BASIS_KEY_NAME         + forCcyName;
    mdmKeys[OswDomModelKey]         = OSWMODEL_KEY_NAME         + domCcyName;
    mdmKeys[OswForModelKey]         = OSWMODEL_KEY_NAME         + forCcyName;
    mdmKeys[FxModelKey]             = FXMODEL_KEY_NAME          + fxName;
    mdmKeys[CorrelMatrixKey]        = CORREL_KEY_NAME           + domForFxName;
    mdmKeys[MrsDomKey]              = MRS_KEY_NAME              + domCcyName;
    mdmKeys[MrsForKey]              = MRS_KEY_NAME              + forCcyName;
    mdmKeys[QFxKey]                 = Q_KEY_NAME                + fxName;
    mdmKeys[QDomKey]                = Q_KEY_NAME                + domCcyName;
    mdmKeys[QForKey]                = Q_KEY_NAME                + forCcyName;

    SetKeys(mdmKeys);

	ComputeProductVectorsFromCurves();

	/// Build an internal cst manager to save
     ARM_CstManagerPtr fxDataManager = ComputeCstManager();

    /// Create the Generic Security paid in domestic=coupon currency
    CreateAndSetDealDescriptionAndTimeIt(GetKeys()[YcBasisDomKey],itsColumnsToPrice,fxDataManager);

    /// To boost time, disable intermediatePayoffs & snapshots computation
    GetGenSecurity()->SetOtherPayoffsFlag(false);
}

	/////////////////////////////////////////////////////////////////
///	Class  : ARM_CCSCalculator
///	Routine: init
///	Returns: void
///	Action : initialize the ccswap calculator
/////////////////////////////////////////////////////////////////
void ARM_CCSCalculator::Init(	const ARM_MarketData_ManagerRep& mktDataManager,
		const ARM_GP_Vector& schedulerDatas,
		const ARM_GP_Vector& truncatorDatas,
		bool markovianDriftSamplerFlag,
		ARM_PRDCCalibType calibType,
		const ARM_GP_Vector& calibDatas)
{
	/// Register input objects but no more than the internal number of access keys
    size_t nbToReg = GetKeys().size();
    for(size_t i(0); i<nbToReg; ++i)
	{
		if(!mktDataManager.TestIfKeyMissing(GetKeys()[i]))
			GetMktDataManager()->RegisterData(GetKeys()[i],mktDataManager.GetData(GetKeys()[i]));
	}
	GetMktDataManager()->SetDetailMode(mktDataManager.GetDetailMode());

	/// Check market datas
    CheckMktDataAndTimeIt();

	itsSchedulerDatas = schedulerDatas;
    itsTruncatorDatas = truncatorDatas;
    itsCalibType = calibType;
    itsCalibDatas = calibDatas;

	/// Create a 2IR+FX model
    CreateAndSetModelAndTimeIt();

	/// Create calibration sets
    CreateAndSetCalibrationAndTimeIt();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CCSCalculator
///	Routine: ComputeProductVectorsFromCurves
///	Returns: void
///	Action : compute vector attributes from ARM_Curve 
/////////////////////////////////////////////////////////////////
void ARM_CCSCalculator::ComputeProductVectorsFromCurves()
{
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();	
	
	/// Convert ARM_Curves into relevant vectors
	/// this will be much easier for computations
	double lag,margin,nominal;

	ARM_GP_Vector* forResetDates = itsForDateStrip->GetResetDates() ;
	for ( size_t i(0); i<itsForSize; i++)
	{
		lag = (*forResetDates)[i]-asOfDate;
		margin = itsForSpreadCv.Interpolate(lag);
		itsvForSpread.push_back(margin );
		nominal = itsForNominalCv.Interpolate(lag);
		itsvForNominal.push_back(nominal);
	}

	ARM_GP_Vector* domResetDates = itsDomDateStrip->GetResetDates() ;
	
	for (i=0; i<itsDomSize; i++)
	{
		lag = (*domResetDates)[i]-asOfDate;
		margin = itsDomSpreadCv.Interpolate(lag);
		itsvDomSpread.push_back(margin );
		nominal = itsDomNominalCv.Interpolate(lag);
		itsvDomNominal.push_back(nominal);

		if (i < itsNbFixFlows)
		{
			itsvCpnIsFixed.push_back(true),
			itsvFixCpn.push_back(0.0);
		}
		else
		{	
			itsvCpnIsFixed.push_back(false);
		}
	}
	
	ARM_GP_Vector* exerciseDates = itsExerciseDateStrip->GetResetDates() ;
	if(!itsFees.empty())
	{
		double fee; 
		for (i=0; i<itsExerSize; i++)
		{
			lag = (*exerciseDates)[i] - asOfDate;
			fee =  itsFees.Interpolate(lag);
			((fee > NON_CALL_FEE ) || (i < itsNbNoCall)) ? itsvIsExerDate.push_back(false) : itsvIsExerDate.push_back(true);
		}
		//itsNbNoCall = CC_NS(std,distance)(itsvIsExerDate.begin(),itsvIsExerDate.find(true));
	}

	/// last: compute indexes to relate exer - funding - cpn		
	ARM_GP_Vector* domStartDates  = itsDomDateStrip->GetFlowStartDates() ;
	ARM_GP_Vector* forStartDates = itsForDateStrip->GetFlowStartDates() ;
	
	i = 0;
	for (size_t k(0); k<itsExerSize; k++)
	{
		while ((*exerciseDates)[k]>(*forStartDates)[i]) i++;
		if((*exerciseDates)[k]>(*forResetDates)[i])
			ARM_THROW( ERR_INVALID_ARGUMENT, "CCS calculator : no call date allowed between fixing dates and start dates");
		itsvForIndex.push_back(i);
	}
	for (k=0, i=0; k<itsExerSize; k++)
	{
		while ((*exerciseDates)[k]>(*domStartDates)[i]) i++;
		if((*exerciseDates)[k]>(*domResetDates)[i])
			ARM_THROW( ERR_INVALID_ARGUMENT, "CCS calculator : no call date allowed between fixing dates and start dates");
		itsvDomIndex.push_back(i);
	}
	itsvForIndex.push_back(forStartDates->size());
	itsvDomIndex.push_back(domStartDates->size());
		
	/* 2 call dates cannot belong to the same period */
	/* this may happen when freq (call) > freq (exotic) or freq (funding) */
	for (k=0; k<itsExerSize; k++)
	{
		if (itsvForIndex[k]>=itsvForIndex[k+1])
			ARM_THROW( ERR_INVALID_ARGUMENT, "CCS calculator : 2 call dates found within the same funding leg period. Freq (call) > Freq (Funding) is forbidden");

		if (itsvDomIndex[k]>=itsvDomIndex[k+1])
			ARM_THROW( ERR_INVALID_ARGUMENT, "CCS calculator : 2 call dates found within the same funding leg period. Freq (call) > Freq (Coupon) is forbidden");
	}

	/* after each call call date, we must have same start dates for exo and funding legs */
	for (k=0; k<itsExerSize; k++) 
	{		
		if (fabs((*forStartDates)[itsvForIndex[k]]-(*domStartDates)[itsvDomIndex[k]])> ARM_GlobalConstant::ARM_SEVENDAYS_LAG)
			ARM_THROW( ERR_INVALID_ARGUMENT, "CCS calculator : call date  defines a swap with mismatching exo and floating start dates (no broken period allowed)");
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CCSCalculator
///	Routine: ComputeCstManager
///	Returns: void
///	Action : Build a constant manager with Fx option date strip &
///          profiles
/////////////////////////////////////////////////////////////////
ARM_CstManagerPtr ARM_CCSCalculator::ComputeCstManager()
{
    ARM_StringVector names(0);
    vector <ARM_GramFctorArg> values(0);

    /// foreign  stripDates
    names.push_back(CCSProfileNamesTable[ForDateStripProfile]);
    values.push_back(ARM_GramFctorArg(itsForDateStrip));

	/// domestic stripDates
    names.push_back(CCSProfileNamesTable[DomDateStripProfile]);
    values.push_back(ARM_GramFctorArg(itsDomDateStrip));
	
	/// foreign margin times notional
    names.push_back(CCSProfileNamesTable[ForMarginTimesNominalProfile]);
    values.push_back(ARM_GramFctorArg(ARM_GP_VectorPtr(new ARM_GP_Vector(itsvForSpread * itsvForNominal) )));

	// domestic margin times notional
    names.push_back(CCSProfileNamesTable[DomMarginTimesNominalProfile]);
    values.push_back(ARM_GramFctorArg(ARM_GP_VectorPtr(new ARM_GP_Vector(itsvDomSpread * itsvDomNominal) )));

	names.push_back(CCSProfileNamesTable[DomNominalProfile]);
	values.push_back(ARM_GramFctorArg(ARM_GP_VectorPtr((ARM_GP_Vector*)itsvDomNominal.Clone() )));

	names.push_back(CCSProfileNamesTable[ForNominalProfile]);
	values.push_back(ARM_GramFctorArg(ARM_GP_VectorPtr((ARM_GP_Vector*)itsvForNominal.Clone() )));

	names.push_back(CCSProfileNamesTable[ForMarginProfile]);
	values.push_back(ARM_GramFctorArg(ARM_GP_VectorPtr((ARM_GP_Vector*)itsvForSpread.Clone() )));

	names.push_back(CCSProfileNamesTable[DomMarginProfile]);
	values.push_back(ARM_GramFctorArg(ARM_GP_VectorPtr((ARM_GP_Vector*)itsvDomSpread.Clone() )));

	return  ARM_CstManagerPtr(new ARM_CstManager(names, values));
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CCSCalculator
///	Routine: CheckData
///	Returns: void
///	Action : check if CCS datas are consistent
/////////////////////////////////////////////////////////////////
void ARM_CCSCalculator::CheckData()
{
    /// Check if input columns to price are present in the deal description
    size_t colNamesSize = sizeof(CCSColNamesTable)/sizeof(CCSColNamesTable[0]);
    for(size_t i=0;i<itsColumnsToPrice.size();++i)
    {
        for(size_t j=0;j<colNamesSize;++j)
            if(itsColumnsToPrice[i] == CCSColNamesTable[j])
                break;
        if(j==itsColumnsToPrice.size())
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : can't price this unknown column name " + itsColumnsToPrice[i]);  
    }
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CCSCalculator
///	Routine: CheckMktData
///	Returns: void
///	Action : check if PRDC market datas are consistent
/////////////////////////////////////////////////////////////////
void ARM_CCSCalculator::CheckMktData()
{
    /// Market datas checking

	ARM_ZeroCurve* domCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcDomKey]));
    if(!domCurve)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Yield Curve for key=" + GetKeys()[YcDomKey] + " is expected in the Market Data Manager");
    string domCcy(domCurve->GetCurrencyUnit()->GetCcyName());

	ARM_ZeroCurve* basisDomCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcBasisDomKey]));
    if(!basisDomCurve)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Yield Curve for key=" + GetKeys()[YcBasisDomKey] + " is expected in the Market Data Manager");
    string basisDomCcy(basisDomCurve->GetCurrencyUnit()->GetCcyName());

    if(domCcy != basisDomCcy)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Domestic(=coupon) Basis Curve currency should consistent with reference curve");

	ARM_ZeroCurve* forCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcForKey]));
    if(!forCurve)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Yield Curve for key=" + GetKeys()[YcForKey] + " is expected in the Market Data Manager");
    string forCcy(forCurve->GetCurrencyUnit()->GetCcyName());

	ARM_ZeroCurve* basisForCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcBasisForKey]));
    if(!basisForCurve)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Yield Curve for key=" + GetKeys()[YcBasisForKey] + " is expected in the Market Data Manager");
    string basisForCcy(basisForCurve->GetCurrencyUnit()->GetCcyName());

    if(forCcy != basisForCcy)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Foreign Basis Curve currency should consistent with reference curve");

	ARM_Forex* forex = dynamic_cast<ARM_Forex*>(GetMktDataManager()->GetData(GetKeys()[ForexKey]));
	if(!forex )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : forex=" + GetKeys()[ForexKey] + " is expected in the Market Data Manager");

    ARM_BSModel* oswDomBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswDomModelKey]) );
    if(!oswDomBSModel)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : domestic swaption B&S model for key=" + GetKeys()[OswDomModelKey] + " is expected in the Market Data Manager");

    ARM_BSModel* oswForBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswForModelKey]) );
    if(!oswForBSModel)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : foreign swaption B&S model for key=" + GetKeys()[OswForModelKey] + " is expected in the Market Data Manager");

	ARM_ModelParam* mrsDomParam = dynamic_cast< ARM_ModelParam* >( GetMktDataManager()->GetData(GetKeys()[MrsDomKey]) );
    if(!mrsDomParam || mrsDomParam->GetType() != ARM_ModelParamType::MeanReversion)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : MRS Param for key=" + GetKeys()[MrsDomKey] + " is expected in the Market Data Manager");

	ARM_ModelParam* mrsForParam = dynamic_cast< ARM_ModelParam* >( GetMktDataManager()->GetData(GetKeys()[MrsForKey]) );
    if(!mrsForParam || mrsForParam->GetType() != ARM_ModelParamType::MeanReversion)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Foreign MRS Param for key=" + GetKeys()[MrsForKey] + " is expected in the Market Data Manager");
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CCSCalculator
///	Routine: ColumnNames
///	Returns: ARM_RowInfo
///	Action : create the column names of the deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_CCSCalculator::ColumnNames() const
{
    size_t colNamesSize = sizeof(CCSColNamesTable)/sizeof(CCSColNamesTable[0]);
    vector< string > colNamesVec(colNamesSize);
    vector< ARM_GP_VALUE_TYPE > colTypeVec(colNamesSize, ARM_STRING); 

    for(size_t i=0;i<colNamesSize; ++i)
        colNamesVec[i] = CCSColNamesTable[i];

    ARM_RowInfo rowInfo(colNamesVec,colTypeVec);

    return rowInfo;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CCSCalculator
///	Routine: CustomDatesStructure
///	Returns: ARM_DateStripCombiner
///	Action : create the list of all event/notice dates of the PRDC.
///			customized Schedule Not implemented yet.
/////////////////////////////////////////////////////////////////
ARM_DateStripCombiner ARM_CCSCalculator::CustomDatesStructure(const ARM_DateStripVector& dateStrips) const
{
	return ARM_DateStripCombiner();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CCSCalculator
///	Routine: DatesStructure
///	Returns: ARM_DateStripCombiner
///	Action : create the list of all event/notice dates of the PRDC.
/////////////////////////////////////////////////////////////////
ARM_DateStripCombiner ARM_CCSCalculator::DatesStructure() const
{
    double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	
	int fwdRule			= K_MOD_FOLLOWING;	// for forward dates
	int intRule			= K_ADJUSTED;			// for interest dates
	int resetTiming		= K_ADVANCE;
	int payTiming		= K_ARREARS;
	
	ARM_INDEX_TYPE indexType = ((ARM_Currency&)itsDomesticCcy).GetVanillaIndexType();
	char* DefaultresetCalendar = itsDomesticCcy.GetResetCalName(indexType);
	CC_NS(std,auto_ptr)<char> holdresetCalendar(DefaultresetCalendar);
	char* DefaultpayCalendar  = itsDomesticCcy.GetPayCalName(indexType);
	CC_NS(std,auto_ptr)<char> holdpayCalendar(DefaultpayCalendar);

	const char* resetCalendar	= itsResetCalendar == "" ? DefaultresetCalendar : itsResetCalendar.c_str();
	const char* payCalendar		= itsPayCalendar == "" ? DefaultpayCalendar: itsPayCalendar.c_str();

	size_t nbPastNoCall;
	ARM_GP_Vector* resetDates = NULL;
	int indexFirstNotice =0;
	size_t size;
	/// 1- exercise schedule
	/// generate schedules only once ...
	if (itsExerciseDateStrip == ARM_DateStripPtr(NULL))
	{
		ARM_DateStrip ExerSched(itsStartDate,itsEndDate,itsExerciseFreq,itsDomDaycount,resetCalendar,fwdRule,intRule,itsStubRule,
			-fabs(itsNoticeGap),itsExerciseFreq, GETDEFAULTVALUE, payCalendar, resetTiming, payTiming);
	
		/// past  no call
		nbPastNoCall=0;
	    resetDates = ExerSched.GetResetDates() ;
		size = ExerSched.size();
		while(nbPastNoCall < size && (*resetDates)[nbPastNoCall] < asOfDate + K_NEW_DOUBLE_TOL)
			++nbPastNoCall;		
		const_cast< ARM_CCSCalculator* >(this)->itsNbNoCall= itsNbNoCall > nbPastNoCall ? itsNbNoCall - nbPastNoCall : 0;

		/// Keep only the futur notification dates
		ExerSched.ResizeAndBuilt(nbPastNoCall,size);

		const_cast< ARM_CCSCalculator* >(this)->itsExerciseDateStrip = ARM_DateStripPtr(new ARM_DateStrip(ExerSched));	
		const_cast< ARM_CCSCalculator* >(this)->itsExerSize = ExerSched.size();		
	}
	double firstNoticeDate = (*itsExerciseDateStrip->GetResetDates())[0];
	
		/// 3- structured schedule
	/// generate schedules only once ...
	if (itsDomDateStrip == ARM_DateStripPtr(NULL))
	{
		/// 4- structured schedule
		ARM_DateStrip domesticSched(itsStartDate,itsEndDate,itsDomFreq,itsDomDaycount,resetCalendar,fwdRule,intRule,
									 itsStubRule,GETDEFAULTVALUE,itsDomFreq, GETDEFAULTVALUE, payCalendar, resetTiming, payTiming);

		nbPastNoCall=0;
		resetDates = domesticSched.GetResetDates() ;
		size = domesticSched.size();
		while(nbPastNoCall < size && (*resetDates)[nbPastNoCall] < asOfDate + K_NEW_DOUBLE_TOL)
			++nbPastNoCall;

		domesticSched.ResizeAndBuilt(nbPastNoCall,size);
		const_cast< ARM_CCSCalculator* >(this)->itsDomDateStrip = ARM_DateStripPtr( new ARM_DateStrip(domesticSched));

		if(fabs(itsFxResetGap) > fabs(itsNoticeGap) + K_NEW_DOUBLE_TOL)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Fx Reset gap must be less than notification gap, please advise");

		/// 4'- forex schedule
		ARM_DateStrip ForexSched(itsStartDate,itsEndDate,itsDomFreq,itsDomDaycount,resetCalendar,fwdRule,intRule,itsStubRule,
								-fabs(itsFxResetGap) ,itsDomFreq, GETDEFAULTVALUE, payCalendar, payTiming, payTiming);

		ForexSched.ResizeAndBuilt(nbPastNoCall,ForexSched.size());
		const_cast< ARM_CCSCalculator* >(this)->itsForexDateStrip = ARM_DateStripPtr( new ARM_DateStrip(ForexSched));
	}
	const_cast< ARM_CCSCalculator* >(this)->itsDomSize = itsDomDateStrip->size();
	double firstCancelDate = (*itsDomDateStrip->GetResetDates())[0];
		
	/// 2- funding schedule

	indexType = ((ARM_Currency&)itsForeignCcy).GetVanillaIndexType();
	char* DefaultforResetCalendar = itsForeignCcy.GetResetCalName(indexType);
	CC_NS(std,auto_ptr)<char> holdforResetCalendar(DefaultforResetCalendar);
	char* DefaultforPayCalendar  = itsDomesticCcy.GetPayCalName(indexType);
	CC_NS(std,auto_ptr)<char> holdForPayCalendar(DefaultforPayCalendar);

	resetCalendar	= itsForResetCalendar == "" ? DefaultforResetCalendar : itsForResetCalendar.c_str();
	payCalendar		= itsForPayCalendar == "" ? DefaultforPayCalendar: itsForPayCalendar.c_str();

	/// generate schedules only once ...
	if (itsForDateStrip == ARM_DateStripPtr(NULL))
	{
		ARM_DateStrip foreignSched(	itsStartDate,itsEndDate,itsForFreq,itsForDaycount,resetCalendar,fwdRule,intRule,
								    itsStubRule,GETDEFAULTVALUE,itsForFreq, GETDEFAULTVALUE, payCalendar, resetTiming,
									payTiming);

		nbPastNoCall=0;
		resetDates = foreignSched.GetResetDates() ;
		size = foreignSched.size();
		while( nbPastNoCall < size  && (*resetDates)[nbPastNoCall] < asOfDate + K_NEW_DOUBLE_TOL)
			++nbPastNoCall;
		while( nbPastNoCall < size  && (*resetDates)[nbPastNoCall] < CC_Min(firstNoticeDate,firstCancelDate)  - ARM_GlobalConstant::ARM_SEVENDAYS_LAG)
			++nbPastNoCall;
		foreignSched.ResizeAndBuilt(nbPastNoCall,size);

		const_cast< ARM_CCSCalculator* >(this)->itsForDateStrip = ARM_DateStripPtr(new ARM_DateStrip(foreignSched));
	}
	
	const_cast< ARM_CCSCalculator* >(this)->itsForSize = itsForDateStrip->size();	

    /// Merge schedules on "ResetDate"
    ARM_DateStripVector SchedVect(1);
    SchedVect[0] = &*itsExerciseDateStrip;

    ARM_DateStripCombiner EventSchedule(SchedVect,"ResetDate");

	return EventSchedule;
}
/////////////////////////////////////////////////////////////////
///	Class  : ARM_CCSCalculator
///	Routine: InitPriceableColumns
///	Returns: void
///	Action : initialise to 0 column to be able to sum
/////////////////////////////////////////////////////////////////
void ARM_CCSCalculator::InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const
{
    string zeroValue("0");

    rowDescVec[FirstCCSSwap] = zeroValue;
    rowTypeVec[FirstCCSSwap] = ARM_DOUBLE;

    rowDescVec[CCSBermuda] = zeroValue;
    rowTypeVec[CCSBermuda] = ARM_DOUBLE;

	rowDescVec[CCSOption] = zeroValue;
    rowTypeVec[CCSOption] = ARM_DOUBLE;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CCSCalculator
///	Routine: MiddleRows
///	Returns: ARM_RowInfo
///	Action : create a row of a deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_CCSCalculator::MiddleRows( size_t eventIdx, const ARM_DateStripCombiner& datesStructure ) const
{
	size_t eventSize = datesStructure.GetDateStrip(0)->GetResetDates()->size();
	size_t descSize = sizeof(CCSColNamesTable)/sizeof(CCSColNamesTable[0]);

	vector< string > rowDescVec(descSize);
	vector< ARM_GP_VALUE_TYPE > rowTypeVec(descSize, ARM_MISSING_TYPE); 

	/// EventDate description
	double eventDate=(*(datesStructure.GetMergeData()))[eventIdx];
	CC_Ostringstream eventDateDesc;
	eventDateDesc << CC_NS(std,fixed) << eventDate;
	rowDescVec[EventDate] = eventDateDesc.str();
	rowTypeVec[EventDate] = ARM_DATE_TYPE;

	/// Set default 0 value for each column to be able to sum it
	ARM_CCSCalculator::InitPriceableColumns(rowDescVec,rowTypeVec);

	/// Get the model names for domestic and forex models
	string domModelName = GetKeys()[YcBasisDomKey];
	string forModelName = GetKeys()[YcBasisForKey];
	string fxModelName = GetKeys()[ForexKey];

	string nextExerIdx("[i+1]");
	bool isLastEvent = (eventIdx+1>=eventSize);

	/// Funding description : in the booster table the spread is already converted
	/// in the domestic ccy then the leg is described as a simple (but fictive)
	/// domestic floating leg

	double date=(*(itsDomDateStrip->GetFlowStartDates()))[itsvDomIndex[eventIdx]];
	CC_Ostringstream domStartDesc;
	domStartDesc << CC_NS(std,fixed) << date;
	rowDescVec[DomStartDate] = domStartDesc.str();
	rowTypeVec[DomStartDate] = ARM_DATE_TYPE;

	date=(*(itsDomDateStrip->GetFlowEndDates()))[itsvDomIndex[eventIdx]];
	CC_Ostringstream domEndDesc;
	domEndDesc << CC_NS(std,fixed) << date;
	rowDescVec[DomEndDate] = domEndDesc.str();
	rowTypeVec[DomEndDate] = ARM_DATE_TYPE;

	date=(*(itsForDateStrip->GetFlowStartDates()))[itsvForIndex[eventIdx]];
	CC_Ostringstream forStartDesc;
	forStartDesc << CC_NS(std,fixed) << date;
	rowDescVec[ForStartDate] = forStartDesc.str();
	rowTypeVec[ForStartDate] = ARM_DATE_TYPE;

	date=(*(itsForDateStrip->GetFlowEndDates()))[itsvForIndex[eventIdx]];
	CC_Ostringstream forEndDesc;
	forEndDesc << CC_NS(std,fixed) << date;
	rowDescVec[ForEndDate] = forEndDesc.str();
	rowTypeVec[ForEndDate] = ARM_DATE_TYPE;	
	
	CC_Ostringstream EndDateDesc;
	EndDateDesc << CC_NS(std,fixed) << itsEndDate.GetJulian();
	rowDescVec[EndDate] = EndDateDesc.str();
	rowTypeVec[EndDate] = ARM_DATE_TYPE;

	CC_Ostringstream forNotionalDesc;
	forNotionalDesc << CC_NS(std,fixed) << itsvForNominal[eventIdx];
	rowDescVec[ForNominal] = forNotionalDesc.str();
	rowTypeVec[ForNominal] = ARM_DOUBLE_TYPE;

	CC_Ostringstream domNotionalDesc;
	domNotionalDesc << CC_NS(std,fixed) << itsvDomNominal[eventIdx];
	rowDescVec[DomNominal] = domNotionalDesc.str();
	rowTypeVec[DomNominal] = ARM_DOUBLE_TYPE;

	/// to convert from long to string
	string forDayCount	= ARM_ArgConvReverse_DayCount.GetString(itsForDaycount);
	string forFreq		= ARM_ArgConvReverse_MatFrequency.GetString(itsForFreq);
	
	/// to convert from long to string
	string domDayCount	= ARM_ArgConvReverse_DayCount.GetString(itsDomDaycount);
	string domFreq		= ARM_ArgConvReverse_MatFrequency.GetString(itsDomFreq);

	//Foreign leg with exchange notional
	string recPay = itsPayRec == K_RCV ? "REC" : "PAY";
	CC_Ostringstream BasisSwapDesc;
	BasisSwapDesc << "BASISSWAP(" << domModelName << ","<< forModelName <<","<< fxModelName <<"," << CCSColNamesTable[ForStartDate] << "[i], "; 
	BasisSwapDesc << CCSColNamesTable[EndDate] << "[i]," << recPay << "," << domFreq << "," << domDayCount << "," ;
	BasisSwapDesc << forFreq << "," << forDayCount << ",FLOTTANT,FLOTTANT,10,2,"<< CCSProfileNamesTable[DomMarginProfile] <<",,";
	BasisSwapDesc << CCSProfileNamesTable[ForMarginProfile] <<",," << CCSProfileNamesTable[DomNominalProfile] <<",";
	BasisSwapDesc << CCSProfileNamesTable[ForNominalProfile]<< ",,,," << itsvDomIndex[eventIdx] << ",," << itsvForIndex[eventIdx] <<")";
	rowDescVec[BasisCCSSwap] = BasisSwapDesc.str();
	rowTypeVec[BasisCCSSwap] = ARM_STRING;
    
	CC_Ostringstream domesticDesc;
	domesticDesc << "SWAP(" << domModelName << "," << CCSColNamesTable[DomStartDate] << "[i], ";
	domesticDesc << CCSColNamesTable[DomEndDate] << "[i], 0, PAY,,,"<< domFreq << ", " << domDayCount << ",";
	domesticDesc << CCSProfileNamesTable[DomMarginProfile] <<","<< CCSProfileNamesTable[DomNominalProfile] <<",,,," ;
	domesticDesc << CCSProfileNamesTable[DomDateStripProfile] << ",," << itsvDomIndex[eventIdx] << ")";
	rowDescVec[DomesticLeg] = domesticDesc.str();
	rowTypeVec[DomesticLeg] = ARM_STRING;

	CC_Ostringstream foreignDesc;
	foreignDesc << "SWAP(" << forModelName << "," << CCSColNamesTable[ForStartDate] << "[i], ";
	foreignDesc << CCSColNamesTable[ForEndDate] << "[i], 0, PAY,,,"<< forFreq << ", " << forDayCount << ",";
	foreignDesc << CCSProfileNamesTable[ForMarginProfile] <<","<< CCSProfileNamesTable[ForNominalProfile] <<",,,," ;
	foreignDesc << CCSProfileNamesTable[ForDateStripProfile] << ",," << itsvForIndex[eventIdx] << ")";
	rowDescVec[ForeignLeg] = foreignDesc.str();
	rowTypeVec[ForeignLeg] = ARM_STRING;

	CC_Ostringstream spotDesc;
	spotDesc << "SPOT(" << fxModelName << ")";
	rowDescVec[Spot] = spotDesc.str();
	rowTypeVec[Spot] = ARM_STRING;

	CC_Ostringstream domNExDesc;
	domNExDesc << CCSColNamesTable[DomNominal] << "[i]*(DF(" << domModelName <<","<< CCSColNamesTable[DomEndDate] << "[i])-";
	domNExDesc << "DF(" << domModelName <<","<< CCSColNamesTable[DomStartDate] << "[i]))";
	rowDescVec[DomExchangeNotional] = domNExDesc.str();
	rowTypeVec[DomExchangeNotional] = ARM_STRING;

	CC_Ostringstream forNExDesc;
	forNExDesc << CCSColNamesTable[ForNominal] << "[i]*(DF(" << forModelName <<","<< CCSColNamesTable[ForEndDate] << "[i])-";
	forNExDesc << "DF(" << forModelName <<","<< CCSColNamesTable[ForStartDate] << "[i]))";
	rowDescVec[ForExchangeNotional] = forNExDesc.str();
	rowTypeVec[ForExchangeNotional] = ARM_STRING;

	CC_Ostringstream foriegnFlowDesc;
	foriegnFlowDesc <<"("<< CCSColNamesTable[ForeignLeg] << "[i]+" << CCSColNamesTable[ForExchangeNotional] << "[i])*";
	foriegnFlowDesc << CCSColNamesTable[Spot]  << "[i]";
	rowDescVec[ForeignFlow] = foriegnFlowDesc.str();
	rowTypeVec[ForeignFlow] = ARM_STRING;

	CC_Ostringstream domesticFlowDesc;
	domesticFlowDesc << CCSColNamesTable[DomesticLeg] << "[i]+" << CCSColNamesTable[DomExchangeNotional] << "[i]";
	rowDescVec[DomesticFlow] = domesticFlowDesc.str();
	rowTypeVec[DomesticFlow] = ARM_STRING;

	CC_Ostringstream ccsFlowDesc;
	if(itsPayRec == K_RCV)
		ccsFlowDesc << CCSColNamesTable[DomesticFlow] << "[i]-" << CCSColNamesTable[ForeignFlow] << "[i]";
	else
		ccsFlowDesc << CCSColNamesTable[ForeignFlow] << "[i]-" << CCSColNamesTable[DomesticFlow] << "[i]";
	rowDescVec[CCSSwapFlow] = ccsFlowDesc.str();
	rowTypeVec[CCSSwapFlow] = ARM_STRING;

	CC_Ostringstream ccSwapDesc;
	if(isLastEvent)
		ccSwapDesc << CCSColNamesTable[CCSSwapFlow] << "[i]";
	else
		ccSwapDesc << CCSColNamesTable[CCSSwapFlow] << "[i]+PV(" << CCSColNamesTable[CCSwap] << "[i+1])";
	rowDescVec[CCSwap] = ccSwapDesc.str();
	rowTypeVec[CCSwap] = ARM_STRING;


	/// CCS option description
	CC_Ostringstream ccsBermudaDesc;
	if(itsvIsExerDate[eventIdx])
	{
		/// Actual exercise at this event date
		ccsBermudaDesc << "EXERCISE(0," << CCSColNamesTable[BasisCCSSwap] << "[i],";
		if(isLastEvent)
			/// No residual option
			ccsBermudaDesc << "0)";
		else
			/// Only one description to get the bermuda price
			ccsBermudaDesc << CCSColNamesTable[CCSBermuda] << nextExerIdx << ")";
	}
	else
	{
		/// No exercise allowed at this event date
		if(isLastEvent)
			/// No residual option
			ccsBermudaDesc << "0";
		else
			/// Just actualise residual option
			ccsBermudaDesc << "PV(" << CCSColNamesTable[CCSBermuda] << nextExerIdx << ")";
	}

	rowDescVec[CCSBermuda] = ccsBermudaDesc.str();
	rowTypeVec[CCSBermuda] = ARM_STRING;

	/// The single CCS option a first line to get the price
	bool isFirstEvent = (eventIdx==itsNbNoCall);
	if(isFirstEvent)
	{
		CC_Ostringstream ccsFirstSwapDesc;
        ccsFirstSwapDesc << CCSColNamesTable[BasisCCSSwap] << "[i]";
        rowDescVec[FirstCCSSwap] = ccsFirstSwapDesc.str();
        rowTypeVec[FirstCCSSwap] = ARM_STRING;

		CC_Ostringstream ccsOptionDesc;
		ccsOptionDesc << CCSColNamesTable[CCSBermuda] << "[i]";
		rowDescVec[CCSOption] = ccsOptionDesc.str();
		rowTypeVec[CCSOption] = ARM_STRING;

	}
	return ARM_RowInfo(rowDescVec,rowTypeVec);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CCSCalculator
///	Routine: CreateAndSetModel
///	Returns: void
///	Action : create the 2IR+FX model
/////////////////////////////////////////////////////////////////
void ARM_CCSCalculator::CreateAndSetModel()
{
    /// Create the 2IR+FX model
    /// Q vol without any curve because will be bootstrapped latter
    ARM_CurveModelParam volParam( ARM_ModelParamType::QVol,0.0,"QVOL");
    ARM_ModelParamVector modelParams(3);
    modelParams[0]  = &volParam;

    /// Multi-assets avec modèles locaux sur le Fx si nécessaire
    int nbModels = 5;

    ARM_StringVector names(nbModels);
	ARM_StringVectorVector depends(nbModels);
    names[ARM_2IRFXModel::DomModel]        = GetKeys()[YcDomKey];
    names[ARM_2IRFXModel::ForModel]        = GetKeys()[YcForKey];
    names[ARM_2IRFXModel::FxModel]         = GetKeys()[ForexKey];	
    names[ARM_2IRFXModel::DomBasisModel]   = GetKeys()[YcBasisDomKey];
    names[ARM_2IRFXModel::ForBasisModel]   = GetKeys()[YcBasisForKey];

	depends[ARM_2IRFXModel::FxModel]       = ARM_StringVector(2);
	depends[ARM_2IRFXModel::FxModel][ARM_2IRFXModel::DomModel]	   = names[ARM_2IRFXModel::DomBasisModel];
	depends[ARM_2IRFXModel::FxModel][ARM_2IRFXModel::ForModel]	   = names[ARM_2IRFXModel::ForBasisModel];
	depends[ARM_2IRFXModel::DomBasisModel] = ARM_StringVector(1,names[ARM_2IRFXModel::DomModel]);
	depends[ARM_2IRFXModel::ForBasisModel] = ARM_StringVector(1,names[ARM_2IRFXModel::ForModel]);
	
    vector< ARM_PricingModelPtr > models(nbModels);
    /// Create the Q1F model for domestic IR market
	modelParams[1] = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(GetKeys()[MrsDomKey]));

	ARM_CurveModelParam QmodelParam = ARM_CurveModelParam(ARM_ModelParamType::QParameter, 0.0,"Q");
	modelParams[2] = &QmodelParam;

	ARM_ZeroCurve* ycDomCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcDomKey]));
    models[ARM_2IRFXModel::DomModel] = ARM_PricingModelPtr( new ARM_QModel1F( CreateClonedPtr(ycDomCurve),ARM_ModelParamsQ1F(modelParams),true) );
	
    /// Create the Q1F model for foreign IR market
	modelParams[1] = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(GetKeys()[MrsForKey]));
	
	ARM_ZeroCurve* ycForCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcForKey]));
    models[ARM_2IRFXModel::ForModel] = ARM_PricingModelPtr( new ARM_QModel1F( CreateClonedPtr(ycForCurve),ARM_ModelParamsQ1F(modelParams),true) );

    /// Create the Q1F model for FX market (with a default MRS=0 because it doesn't matter !)
	ARM_CurveModelParam mrsFxParam( ARM_ModelParamType::MeanReversion,0.0,"FXMRS");
    modelParams[1]   = &mrsFxParam;
    if(GetMktDataManager()->TestIfKeyMissing(GetKeys()[QFxKey]))
    {
        /// Degenerate to pure LN
		QmodelParam = ARM_CurveModelParam(ARM_ModelParamType::QParameter, 0.0,"Q");
        modelParams[2] =  &QmodelParam;
    }
    else
	    modelParams[2] = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(GetKeys()[QFxKey]));

	ARM_ZeroCurve* ycBasisDomCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcBasisDomKey]));
	ARM_ZeroCurve* ycBasisForCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcBasisForKey]));
    ARM_GP_Matrix* correlMatrix = dynamic_cast<ARM_GP_Matrix*>(GetMktDataManager()->GetData(GetKeys()[CorrelMatrixKey]));

    ARM_Forex* forex = dynamic_cast<ARM_Forex*>(GetMktDataManager()->GetData(GetKeys()[ForexKey]));
    models[ARM_2IRFXModel::FxModel] = ARM_PricingModelPtr( ARM_EqFx_ModelFactory.Instance()->CreateModel( CreateClonedPtr(ycBasisDomCurve), 
									modelParams, forex->GetMarketPrice(), 
									CreateClonedPtr(ycBasisForCurve),
									*correlMatrix,
									ARM_EqFx_ModelFactoryImp::Q1F_Model ) );

    /// Create both domestic & foreign forward margin models
    models[ARM_2IRFXModel::DomBasisModel] = ARM_PricingModelPtr( new ARM_ForwardMarginBasis(CreateClonedPtr(ycBasisDomCurve)) );
    models[ARM_2IRFXModel::ForBasisModel] = ARM_PricingModelPtr( new ARM_ForwardMarginBasis(CreateClonedPtr(ycBasisForCurve)) );

	/// Create a modelnamemap
	ARM_ModelNameMap modelMap( names, models, depends );

	/// Finally, We create and set the brand new 2IR+FX Ferrari model !
    ARM_PricingModelPtr hybridModel = ARM_PricingModelPtr( new ARM_2IRFXModel(modelMap,*correlMatrix) );

    /// Create a 3 factors version of the tree ND then set it
    /// The tree is paramerised to replicate the Kernel version
    int schedulerType=ARM_SchedulerBase::MultiRegime;
    int samplerType = itsMarkovianDriftSamplerFlag ? ARM_SamplerBase::MarkovianDrift: ARM_SamplerBase::DriftedMeanReverting;
    int truncatorType=ARM_TruncatorBase::ArrowDebreu;
    int reconnectorType=ARM_ReconnectorBase::Mean;
    int smootherType=ARM_SmootherBase::Linear;
    ARM_TreeBase* tree = ARM_TreeFactory.Instance()->CreateTreeND(3,schedulerType,itsSchedulerDatas,
						samplerType,ARM_GP_Vector(0),truncatorType,itsTruncatorDatas,false,reconnectorType,smootherType);
    hybridModel->SetNumMethod( ARM_NumMethodPtr( tree ) );


    /// Create and set the cash numeraire
    ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::Cash ) );
    hybridModel->SetNumeraire(numeraire);

	/// Set the model
	SetPricingModel(hybridModel);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CCSCalculator
///	Routine: CreateAndSetCalibration
///	Returns: void
///	Action : create the calibration
/////////////////////////////////////////////////////////////////
void ARM_CCSCalculator::CreateAndSetCalibration()
{
    ARM_ModelParamVector emptyCalibParam;
    ARM_CalibMethod* noLinkedMethod=NULL;
    ARM_CalibMethod* noPreviousMethod=NULL;
    bool isCalibShared=false;

    bool isFreezeWeights = false;   // reinit non null weights in the portfolio
    bool isInitVol = true;          // init vol param for calibration (bounds & init guess)

    /// Diagonal domestic swaptions
    /// ---------------------------
    /// Create standard diagonal swaption portfolios (domestic & foreign sigmas calibration)
    pair< ARM_StdPortfolioPtr, ARM_StdPortfolioPtr > diagonalSwaptionPF(CreateDiagonalSwaption());


    /// Build an empty calib method for latter domestic volatility bootstrapping
    ARM_Currency* domCcy = GetCurrencyUnit();
    ARM_CalibMethod* domCalib = CreateIRCalibration(diagonalSwaptionPF.first,domCcy,ARM_2IRFXModel::DomModel);
    /// Compute domestic target prices
    ComputeIROptionPrices(domCalib,OswDomModelKey,isFreezeWeights,isInitVol);

    /// Diagonal foreign swaptions
    /// --------------------------
    /// Build an empty calib method for latter foreign volatility bootstrapping
	ARM_Currency* forCcy = static_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcForKey]))->GetCurrencyUnit();
    ARM_CalibMethod* forCalib = CreateIRCalibration(diagonalSwaptionPF.second,forCcy, ARM_2IRFXModel::ForModel);
    /// Compute foreign target prices
    ComputeIROptionPrices(forCalib,OswForModelKey,isFreezeWeights,isInitVol);

    /// ATM Forex options
    /// -----------------
    /// Create a Fx option portfolio and compute Fx option target prices
    ARM_StdPortfolioPtr fxOptionPF      = CreateFxOption();
    ARM_CalibMethod* fxCalib         = CreateFxCalibration(fxOptionPF,ARM_2IRFXModel::FxModel);
    ComputeFxOptionPrices(fxCalib,isFreezeWeights,isInitVol);

    /// Set next methods to handle calibration with a single method
    domCalib->SetNextMethod(forCalib);
    forCalib->SetNextMethod(fxCalib);

    SetCalibMethod(ARM_CalibMethodPtr(domCalib));
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CCSCalculator
///	Routine: CreateDiagonalSwaption
///	Returns: a portfolio
///	Action : create the list of diagonal swaptions for domestic market
/////////////////////////////////////////////////////////////////
pair< ARM_StdPortfolioPtr, ARM_StdPortfolioPtr > ARM_CCSCalculator::CreateDiagonalSwaption()
{
    double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

    /// Get market swaption standard expiries for domestic & foreign volatilities
    size_t i,nbExp;
    ARM_BSModel* oswBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswDomModelKey]) );
    ARM_VolCurve* oswBSVol = oswBSModel->GetVolatility();
    nbExp = oswBSVol->GetExpiryTerms()->GetSize();
    ARM_GP_Vector domStdExp(nbExp+1);
    domStdExp[0] = asOfDate;
    for(i=0;i<nbExp;++i)
        domStdExp[i+1] = asOfDate + K_YEAR_LEN * (*(oswBSVol->GetExpiryTerms()))[i];

    oswBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswForModelKey]) );
    oswBSVol = oswBSModel->GetVolatility();
    nbExp = oswBSVol->GetExpiryTerms()->GetSize();
    ARM_GP_Vector forStdExp(nbExp+1);
    forStdExp[0] = asOfDate;
    for(i=0;i<nbExp;++i)
        forStdExp[i+1] = asOfDate + K_YEAR_LEN * (*(oswBSVol->GetExpiryTerms()))[i];

    ARM_Swaption* swaption;

    list < ARM_Security* > domSwaptionList;
    list < ARM_Security* > forSwaptionList;


    const ARM_DealDescription& dealDesc = GetGenSecurity()->GetDealDescription();
    size_t nbEvents=dealDesc.GetRowsNb(); // + 1 because of column names
    size_t nbFutureExer=itsvIsExerDate.size();
    if(nbFutureExer+1 != nbEvents)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : mismatch between exercise flags & notice dates in the deal description");

    /// Standard index of the domestic & foreign currencies
    ARM_Currency* domCcy = GetCurrencyUnit();
    int domSpotDays = domCcy->GetSpotDays();
    ARM_INDEX_TYPE domLiborIndex = domCcy->GetVanillaIndexType();
    char* domResetCal = domCcy->GetResetCalName(domLiborIndex);

	ARM_Currency* forCcy = static_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcForKey]))->GetCurrencyUnit();
    int forSpotDays = forCcy->GetSpotDays();
    ARM_INDEX_TYPE forLiborIndex = forCcy->GetVanillaIndexType();
    char* forResetCal = forCcy->GetResetCalName(forLiborIndex);

    /// Diagonal swaptions are based on funding start/end dates but with the right currency
    ARM_Date swapStartDate,expiryDate;
    bool isDomSwap,isForSwap;
    ARM_Swap domSwap,forSwap;

    /// The 1st notice date is always added then get the very last expiry
    double lastNotice,lastInsertedNotice;
    isDomSwap=false;
    size_t firstNoticeIdx=1;
    for(size_t eventIdx=1;eventIdx<nbEvents;++eventIdx)
    {
        if(itsvIsExerDate[eventIdx-1])
        {
            expiryDate = ARM_Date( atof(dealDesc.GetElem(eventIdx,EventDate).c_str()) );
            if(!isDomSwap)
            {
                swapStartDate   = expiryDate;
                swapStartDate.GapBusinessDay(domSpotDays,domResetCal);

                domSwap = ARM_Swap(swapStartDate,itsEndDate,domLiborIndex,0.0,K_MARKET_RATE,K_RCV,
                                 K_DEF_FREQ,K_DEF_FREQ,domCcy);
                swaption = new ARM_Swaption(&domSwap,K_RCV,K_EUROPEAN,K_MARKET_RATE,expiryDate);
                domSwaptionList.push_back(static_cast< ARM_Security* >(swaption));
                isDomSwap=true;

                swapStartDate   = expiryDate;
                swapStartDate.GapBusinessDay(forSpotDays,forResetCal);

                forSwap = ARM_Swap(swapStartDate,itsEndDate,forLiborIndex,0.0,K_MARKET_RATE,K_RCV,
                                 K_DEF_FREQ,K_DEF_FREQ,forCcy);
                swaption = new ARM_Swaption(&forSwap,K_RCV,K_EUROPEAN,K_MARKET_RATE,expiryDate);
                forSwaptionList.push_back(static_cast< ARM_Security* >(swaption));

                firstNoticeIdx = eventIdx;
                lastInsertedNotice = expiryDate.GetJulian();
            }
            lastNotice = expiryDate.GetJulian();
        }
    }
    if(!isDomSwap)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " creation of calibration diagonal swaption not possible");


    /// Between each standard market expiries, only one swaption is created corresponding to
    /// the earlier notice date
    /// To be consistent with Tree3F code (may be it is a bug !), the very last notice date
    /// is excluded
    size_t domExpIdx=0,forExpIdx=0;
    for(eventIdx = firstNoticeIdx+1;eventIdx<nbEvents;++eventIdx)
    {
        expiryDate = ARM_Date( atof(dealDesc.GetElem(eventIdx,EventDate).c_str()) );

        if(itsvIsExerDate[eventIdx-1] && expiryDate.GetJulian() < lastNotice)  /// to be consistent to Tree3F code !!
        {
            isDomSwap=false;
            while(!isDomSwap && domExpIdx < domStdExp.size() && domStdExp[domExpIdx] + K_NEW_DOUBLE_TOL < expiryDate.GetJulian())
            {
                if( ( domExpIdx+1 == domStdExp.size() ||
                      (domExpIdx+1 < domStdExp.size() && expiryDate.GetJulian() < domStdExp[domExpIdx+1] + K_NEW_DOUBLE_TOL) )
                    && lastInsertedNotice + K_NEW_DOUBLE_TOL < domStdExp[domExpIdx] )
                {
                    /// Create a domestic underlying swap
                    isDomSwap=true;
                    swapStartDate   = expiryDate;
                    swapStartDate.GapBusinessDay(domSpotDays,domResetCal);
                    domSwap = ARM_Swap(swapStartDate,itsEndDate,domLiborIndex,0.0,K_MARKET_RATE,K_RCV,
                                     K_DEF_FREQ,K_DEF_FREQ,domCcy);
                }

                ++domExpIdx;
            }

            isForSwap=false;
            while(!isForSwap && forExpIdx < forStdExp.size() && forStdExp[forExpIdx] + K_NEW_DOUBLE_TOL < expiryDate.GetJulian())
            {
                if( ( forExpIdx+1 == forStdExp.size() ||
                      (forExpIdx+1 < forStdExp.size() && expiryDate.GetJulian() < forStdExp[forExpIdx+1] + K_NEW_DOUBLE_TOL) )
                    && lastInsertedNotice + K_NEW_DOUBLE_TOL < forStdExp[forExpIdx] )
                {
                    /// Create a Forestic underlying swap
                    isForSwap=true;
                    swapStartDate   = expiryDate;
                    swapStartDate.GapBusinessDay(forSpotDays,forResetCal);
                    forSwap = ARM_Swap(swapStartDate,itsEndDate,forLiborIndex,0.0,K_MARKET_RATE,K_RCV,
                                     K_DEF_FREQ,K_DEF_FREQ,forCcy);
                }

                ++forExpIdx;
            }


            if(isDomSwap)
            {
                /// Domestic swaption
                swaption = new ARM_Swaption(&domSwap,K_RCV,K_EUROPEAN,K_MARKET_RATE,expiryDate);
                domSwaptionList.push_back(static_cast< ARM_Security* >(swaption));
                lastInsertedNotice = expiryDate.GetJulian();
            }

            /// Foreign swaption
            if(isForSwap)
            {
                swaption = new ARM_Swaption(&forSwap,K_RCV,K_EUROPEAN,K_MARKET_RATE,expiryDate);
                forSwaptionList.push_back(static_cast< ARM_Security* >(swaption));
                lastInsertedNotice = expiryDate.GetJulian();
            }
        }
    }

    ARM_StdPortfolio* domPort = new ARM_StdPortfolio(domSwaptionList);
    for(i=0;i<domPort->size();++i)
    {
        domPort->SetWeight(OSW_DEFAULT_WEIGHT,i);
        domPort->SetPrice((i+1)*OSW_DEFAULT_PRICE,i);
    }
    ARM_StdPortfolio* forPort = new ARM_StdPortfolio(forSwaptionList);
    for(i=0;i<forPort->size();++i)
    {
        forPort->SetWeight(OSW_DEFAULT_WEIGHT,i);
        forPort->SetPrice((i+1)*OSW_DEFAULT_PRICE,i);
    }

    /// Don't forget to free char* calendars !
    delete domResetCal;
    delete forResetCal;


    return pair< ARM_StdPortfolioPtr, ARM_StdPortfolioPtr >(ARM_StdPortfolioPtr(domPort),ARM_StdPortfolioPtr(forPort));
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CCSCalculator
///	Routine: CreateIRCalibration
///	Returns: void
///	Action : create the calibration for IR parameters (vol & MRS)
/////////////////////////////////////////////////////////////////
ARM_CalibMethod* ARM_CCSCalculator::CreateIRCalibration(const ARM_StdPortfolioPtr& diagonalSwaptionPF,
            ARM_Currency* ccy, int modelIdx)
{
    ARM_ModelParamVector emptyCalibParam;
    ARM_CalibMethod* noLinkedMethod=NULL;
    ARM_CalibMethod* noPreviousMethod=NULL;
    bool isCalibShared=false;


    /// Build an empty bootstrap calib method (filled latter for domestic volatility calibration)
    ARM_CalibMethod* volCalib = new ARM_CalibMethod(diagonalSwaptionPF,emptyCalibParam,
                                    ARM_CalibMethodType::Bootstrap1D,ARM_MAX_ITER,
                                    ARM_CalibrationTarget::PriceTarget,
                                    noLinkedMethod,noPreviousMethod,
                                    isCalibShared,modelIdx);

    return volCalib;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CCSCalculator
///	Routine: ComputeIROptionPrice
///	Returns: nothing
///	Action : compute market target prices for calibration purpose
/////////////////////////////////////////////////////////////////
void ARM_CCSCalculator::ComputeIROptionPrices(ARM_CalibMethod* calibMethod, 
	          mdmKeysAlias oswModelIdx, 
			  bool isFreezeWeights, 
			  bool isInitVolParam)
{
    ARM_BSModel* oswBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[oswModelIdx]) );
	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();

    /// Restore calibration portfolios
    ARM_CalibMethod* volCalibMethod;
    ARM_StdPortfolioPtr oswPortfolio;

    double price,vega,nominal,weight;
    ARM_Swaption* swaption;
    size_t i;
    bool isNotSelected;

    volCalibMethod = calibMethod;


    /// Update datas for q volatility calibration
    /// -----------------------------------------
    oswPortfolio = volCalibMethod->GetPortfolio();
    size_t nbOSW = oswPortfolio->GetSize();
    size_t volParamSize=volCalibMethod->GetCalibParams().size();
    bool isInitVol = isInitVolParam || volParamSize == 0;

    ARM_GP_Vector initTimes(nbOSW);
    ARM_GP_Vector initVols(nbOSW);
    double optMat,swapMat,volATM,strike,swapRate;

    for(i=0;i<nbOSW;++i)
    {
        swaption=static_cast< ARM_Swaption* >(oswPortfolio->GetAsset(i));
	    swaption->SetModel(oswBSModel);
        price=swaption->ComputePrice();
        vega=swaption->ComputeSensitivity(K_VEGA);

        nominal = swaption->GetAmount()->CptReferenceValue(0.0); // Cst RefValue
        weight  = (*(oswPortfolio->GetWeights()))[i];

        isNotSelected = vega < IR_VEGA_MIN_TO_SELECT * nominal && !isFreezeWeights;

        oswPortfolio->SetWeight(isNotSelected ? 0.0 : weight,i);
        oswPortfolio->SetPrecision(0.001*vega,i);
        oswPortfolio->SetPrice(price,i);

        if(isInitVol)
        {
            /// Vol initialisation
            optMat  = (swaption->GetExpiryDate().GetJulian()-asOfDate.GetJulian())/K_YEAR_LEN;
            swapMat = (swaption->GetEndDate().GetJulian() - swaption->GetStartDate().GetJulian())/K_YEAR_LEN;
            strike  = swaption->GetStrike();
            volATM  = oswBSModel->ComputeVol(optMat,swapMat,strike,strike)/100.0;

            initTimes[i]    = swaption->GetExpiryDate().GetJulian()-asOfDate.GetJulian();
                swapRate = swaption->CptMarketSwapRate();
                initVols[i] = volATM * swapRate;
        }
    }

    if(isInitVol)
    {
        /// Replace the sigma param with new initialisations
        ARM_GP_Vector volLowerBound(nbOSW,HWVOL_LOWER_BOUND);
        ARM_GP_Vector volUpperBound(nbOSW,HWVOL_UPPER_BOUND);
        ARM_CurveModelParam* vol = new ARM_CurveModelParam(ARM_ModelParamType::QVol,&initVols,&initTimes,
            "QVOL","STEPUPRIGHT",&volLowerBound,&volUpperBound);
        if(volParamSize==0)
            volCalibMethod->GetCalibParams().push_back(vol);
        else if(volParamSize==1)
        {
            delete volCalibMethod->GetCalibParam(0);
            (volCalibMethod->GetCalibParams())[0] = vol;
        }
        else
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : diagonal swaption calibrate only volatilities");
    }
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CCSCalculator
///	Routine: CreateFxOption
///	Returns: a portfolio
///	Action : create the list of ATM FX options
/////////////////////////////////////////////////////////////////
ARM_StdPortfolioPtr ARM_CCSCalculator::CreateFxOption()
{
    double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

	/// To support new FX GP models !
    ARM_BSModel* fxBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[FxModelKey]) );
    
    ARM_Option* fxOption;
    list < ARM_Security* > fxOptionList;

    const ARM_DealDescription& dealDesc = GetGenSecurity()->GetDealDescription();
    size_t i,nbEvents=dealDesc.GetRowsNb(); // + 1 because of column names
    size_t nbFutureExer=itsvIsExerDate.size();
    if(nbFutureExer+1 != nbEvents)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : mismatch between exercise flags & notice dates in the deal description");

    /// Restore forex
	ARM_Forex* forex = static_cast< ARM_Forex* >( GetMktDataManager()->GetData(GetKeys()[ForexKey]) );
	
	ARM_GP_Vector* exerResetDates = itsExerciseDateStrip->GetResetDates();
    if(itsCalibType == ARM_PRCSCalibTypes::ATMCalib)
    {
        /// Portfolio = market standard FX options (ATM)

		/// Compute standard vol expiries
		ARM_VolCurve* fxBSVol = fxBSModel->GetVolatility();
		size_t nbExp = fxBSVol->GetExpiryTerms()->GetSize();
		ARM_GP_Vector fxStdExp(nbExp);
		for(i=0;i<nbExp;++i)
			fxStdExp[i] = asOfDate + K_YEAR_LEN * (*(fxBSVol->GetExpiryTerms()))[i];

        /// Check if the 1st notice date may replace the 1st standard fx option expiry
        ARM_Date firstNoticeDate;
        size_t fxExpIdx=0;
        for(size_t exerIdx=0; exerIdx<itsExerSize; ++exerIdx)
        {
            if(itsvIsExerDate[exerIdx] && !itsvCpnIsFixed [exerIdx])
            {
                firstNoticeDate = ARM_Date( (*itsExerciseDateStrip->GetResetDates())[exerIdx]);
                if(firstNoticeDate.GetJulian() >= fxStdExp[fxExpIdx]-K_NEW_DOUBLE_TOL)
                {
                    /// Replace the 1st standard expiry by the 1st notice
                    fxOption    = new ARM_Option(forex,firstNoticeDate,K_MARKET_RATE,K_CALL,K_EUROPEAN);
                    fxOptionList.push_back(static_cast< ARM_Security* >(fxOption));
                    ++fxExpIdx;
                    break;
                }
            }
        }

        /// Create a fx option for each following standard expiry greater than the 1st notice
        /// but lower or equal to 30y
        double maxFxOptionExpiry = asOfDate + 11111; /// 11111/365 is > 30y!!
        for(;fxExpIdx < nbExp && fxStdExp[fxExpIdx] < maxFxOptionExpiry;++fxExpIdx)
        {
            if(fxStdExp[fxExpIdx] > firstNoticeDate.GetJulian() + K_NEW_DOUBLE_TOL)
            {
                fxOption    = new ARM_Option(forex,ARM_Date(fxStdExp[fxExpIdx]),K_MARKET_RATE,K_CALL,K_EUROPEAN);
                fxOptionList.push_back(static_cast< ARM_Security* >(fxOption));
            }
        }
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
///	Class  : ARM_CCSCalculator
///	Routine: CreateFxCalibration
///	Returns: void
///	Action : create the calibration for IR parameters (vol & MRS)
/////////////////////////////////////////////////////////////////
ARM_CalibMethod* ARM_CCSCalculator::CreateFxCalibration(const ARM_StdPortfolioPtr fxOptionPF, int modelIdx)
{
    ARM_ModelParamVector emptyCalibParam;
    ARM_CalibMethod* noLinkedMethod=NULL;
    ARM_CalibMethod* noPreviousMethod=NULL;
    bool isCalibShared=false;


    /// Build an empty bootstrap calib method (filled latter for fx volatility calibration)
    ARM_CalibMethod* fxCalib = new ARM_CalibMethod(fxOptionPF,emptyCalibParam,
                                    ARM_CalibMethodType::Bootstrap1D,ARM_MAX_ITER,
                                    ARM_CalibrationTarget::PriceTarget,
                                    noLinkedMethod,noPreviousMethod,
                                    isCalibShared,modelIdx);

    return fxCalib;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CCSCalculator
///	Routine: ComputeFxOptionPrice
///	Returns: nothing
///	Action : compute market target prices of the FX option
/////////////////////////////////////////////////////////////////
void ARM_CCSCalculator::ComputeFxOptionPrices(ARM_CalibMethod* calibMethod, bool isFreezeWeights, bool isInitVolParam)
{
	/// To support new FX GP models !
    ARM_BSModel* fxBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[FxModelKey]) );
	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();

    ARM_StdPortfolioPtr fxPortfolio = calibMethod->GetPortfolio();
    size_t nbFx = fxPortfolio->GetSize();
    size_t volParamSize=calibMethod->GetCalibParams().size();
    bool isInitVol = isInitVolParam || volParamSize == 0;

    ARM_GP_Vector initTimes(nbFx);
    ARM_GP_Vector initVols(nbFx);

    ARM_Option* fxOption;
    double optTime,vol,strike,price,vega=0.0,nominal,weight,fwd, atmVol;
    bool isNotSelected;
    size_t strikeIdx=0;
    for(size_t i(0);i<nbFx;++i)
    {
        fxOption = static_cast< ARM_Option* >(fxPortfolio->GetAsset(i));
	    fxOption->SetModel(fxBSModel);
        price=fxOption->ComputePrice();

        switch(itsCalibType)
        {
        case ARM_PRCSCalibTypes::ATSFxMixedCalib : /// Choose ATM/ATS strike with minimum volatility
            vol = fxOption->GetCalcVol();
            strike = fxOption->GetStrike();
            fxOption->SetStrike(K_MARKET_RATE);
            price=fxOption->ComputePrice();
            atmVol = fxOption->GetCalcVol();
            if(vol < atmVol)
            {
                fxOption->SetStrike(strike);
                price=fxOption->ComputePrice();
            }
            break;

        case ARM_PRCSCalibTypes::ATSFxMoneynessCalib :
            if(itsCalibType == ARM_PRCSCalibTypes::ATSFxMoneynessCalib)                 
				///strike = moneyness * fwdFx
                fxOption->SetStrike(itsCalibDatas[strikeIdx] *fxOption->GetCalcFwd());
			else if(itsCalibType == ARM_PRCSCalibTypes::HybridBasketCalib)
				 /// First calibration is ATSFxMoneyness=100% like
                fxOption->SetStrike(fxOption->GetCalcFwd());
            else
                /// strike = shift + fwdFx
                fxOption->SetStrike(itsCalibDatas[strikeIdx] + fxOption->GetCalcFwd());

            price=fxOption->ComputePrice();
            if(strikeIdx+1 < itsCalibDatas.size())
                ++strikeIdx;
            break;

        case ARM_PRCSCalibTypes::ATSFxMinVolCalib : /// Choose ATS with mimimum volatility
            if(itsCalibDatas.size()==0 || i <= itsCalibDatas[0])
            {
                double minMoneyness = FX_MIN_MONEYNESS;
                double maxMoneyness = FX_MAX_MONEYNESS;
                double stepMoneyness= FX_STEP_MONEYNESS;
                if(itsCalibDatas.size()>1)
                    minMoneyness = itsCalibDatas[1];
                if(itsCalibDatas.size()>2)
                    maxMoneyness = itsCalibDatas[2];
                if(itsCalibDatas.size()>3)
                    stepMoneyness = itsCalibDatas[3];

                double moneyness = minMoneyness;
                double minVolStrike = fxOption->GetStrike();
                double minVol = fxOption->GetCalcVol();
                while (moneyness <= maxMoneyness)
                {
                    fxOption->SetStrike(moneyness * fxOption->GetCalcFwd());
                    price=fxOption->ComputePrice();
                    if(fxOption->GetCalcVol() < minVol)
                    {
                        minVolStrike = fxOption->GetStrike();
                        minVol = fxOption->GetCalcVol();
                    }
                    moneyness += stepMoneyness;
                }
                fxOption->SetStrike(minVolStrike);
                price=fxOption->ComputePrice();
            }
            break;
        }

        if(!isFreezeWeights)
            vega=fxOption->ComputeSensitivity(K_VEGA);

        nominal = fxOption->GetAmount()->CptReferenceValue(0.0); // Cst RefValue
        weight  = (*(fxPortfolio->GetWeights()))[i];

        isNotSelected = vega < FX_VEGA_MIN_TO_SELECT * nominal && !isFreezeWeights;

        fxPortfolio->SetWeight(isNotSelected ? 0.0 : weight,i);
        fxPortfolio->SetPrecision(0.001*vega,i);
        fxPortfolio->SetPrice(price,i);

        if(isInitVol)
        {
            /// Vol initialisation (may be too high because forward FX vol is used)
            optTime = fxOption->GetExpiryDate().GetJulian()-asOfDate.GetJulian();
            strike  = fxOption->GetStrike();
            fwd = fxOption->GetCalcFwd();
            vol = fxOption->GetCalcVol()/100.0;

            initTimes[i] = optTime;
            initVols[i]  = vol;
        }
    }

    if(isInitVol)
    {
        /// Replace the sigma param with new initialisations
        ARM_GP_Vector volLowerBound(nbFx,QVOL_LOWER_BOUND);
        ARM_GP_Vector volUpperBound(nbFx,QVOL_UPPER_BOUND);
        ARM_CurveModelParam* vol = new ARM_CurveModelParam(ARM_ModelParamType::QVol,&initVols,&initTimes,
            "QVOL","STEPUPRIGHT",&volLowerBound,&volUpperBound);
        if(volParamSize==0)
            calibMethod->GetCalibParams().push_back(vol);
        else if(volParamSize==1)
        {
            delete calibMethod->GetCalibParam(0);
            (calibMethod->GetCalibParams())[0] = vol;
        }
        else
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : fx options calibrate only volatilities at the moment");
    }
}

////////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCCalculator
///	Routine: ComputePricingData
///	Returns: 
///	Action : pricing function called from the addin GetPricingData()
////////////////////////////////////////////////////////////////////
void ARM_CCSCalculator::ComputePricingData() const
{
	if (!itsHasBeenComputed)
		const_cast<ARM_CCSCalculator*>(this)->PriceAndTimeIt();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CCSCalculator
///	Routine: Calibraton
///	Returns: a double
///	Action : to calibrate the 2IRFxModel. At beginning, we run an auto-calibration 
///          by bootstrapping of the volatility on diagonal swaptions
///          on both domestic & foreign market. Secondly, the spot FX volatility
////         is boostrapped on FX options portfolio. finnaly, we calibrat
///          the local models 
/////////////////////////////////////////////////////////////////
void ARM_CCSCalculator::Calibrate()
{ 
	ARM_2IRFXModel* hybridModel = static_cast< ARM_2IRFXModel* > (&*GetPricingModel());

    /// Calibrate stochatic models and check errors
    GetCalibMethod()->Calibrate(hybridModel);

	/// to check Calibration to avoid any variance squeeze
    CheckCalibErrors();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CCSCalculator
///	Routine: Price
///	Returns: a double
///	Action : price the PRDC option. An auto-calibration is done before
///          by bootstrapping the volatility on diagonal swaptions
///          on both domestic & foreign market then the spot FX volatility
////         is boostrapped on FX options portfolio
/////////////////////////////////////////////////////////////////
double ARM_CCSCalculator::Price()
{
    /// Calibrate
    CalibrateAndTimeIt();
	
	ARM_2IRFXModel* hybridModel = static_cast< ARM_2IRFXModel* > (&*GetPricingModel());

    ARM_PricingModelPtr initialModel;

	double firstColumnPrice = 0.0;
	if (itsColumnsToPrice.size() > 0)
	{
		if(itsMarkovianDriftSamplerFlag)
			/// Save calibrated model because of interpolation over tree schedule
			initialModel = ARM_PricingModelPtr( static_cast< ARM_PricingModel* >(hybridModel->Clone()) );

		ARM_GenPricer* genPricer = new ARM_GenPricer( &*GetGenSecurity(),hybridModel);

		ARM_AutoCleaner<ARM_GenPricer> HoldGP(genPricer );

		firstColumnPrice = genPricer->Price() ;

		string columnName;
		for(size_t i=0;i<itsColumnsToPrice.size();++i)
		{
			columnName = itsColumnsToPrice[i];
			GetPricingData()[columnName] = genPricer->GetPricerInfo()->GetContents(columnName).GetData("Price").GetDouble();
		}
		
		/// Restore calibrated model
		if(itsMarkovianDriftSamplerFlag)
			SetPricingModel(initialModel);
	}
	itsHasBeenComputed = true;

    return firstColumnPrice;       
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CCSCalculator
///	Routine: UpdateCalibration
///	Returns: void
///	Action : update the calibration datas w.r.t. MarketDataManager
///          The context is an hedge ratios computation
/////////////////////////////////////////////////////////////////
void ARM_CCSCalculator::UpdateCalibration(bool isUpdateStrike)
{
	bool isInitVol = true;
	bool isFreezeWeights = false;
    /// Compute domestic target prices
    ComputeIROptionPrices(&*GetCalibMethod(),OswDomModelKey,isFreezeWeights,isInitVol);

	/// Compute foreign target prices
    ComputeIROptionPrices(GetCalibMethod()->GetNextMethod(),OswForModelKey,isFreezeWeights,isInitVol);
    ComputeFxOptionPrices(GetFxCalib(),isFreezeWeights,isInitVol);

	itsHasBeenComputed = false;
}

////////////////////////////////////////////////////
///	Class   : ARM_CCSCalculator
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_CCSCalculator::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream prdcData;

    /// GenCalculator general datas dump
	prdcData << "\n\nCommon Calculator Part\n" << ARM_GenCalculator::toString(indent,nextIndent) << "\n\n";

    return prdcData.str();
}

////////////////////////////////////////////////////
///	Class   : ARM_CCSCalculator
///	Routines: CheckCalibErrors
///	Returns :
///	Action  : Check calibration errors and throw an
///           exception if errors are too large
////////////////////////////////////////////////////
void ARM_CCSCalculator::CheckCalibErrors()
{
    /// Recall that Dom -> For -> Fx with "next" links

    size_t i,nbErrors,calibIdx=0;
    double err;
    ARM_CalibMethod* calib = &*(GetCalibMethod());
    ARM_ModelFitterPtr modelFitter;


    while( calib != NULL &&
           (modelFitter=calib->GetModelFitter()) != ARM_ModelFitterPtr(NULL) )
    {
        modelFitter->SetUpError();

        nbErrors = modelFitter->GetError()->rows();
        for(i=0; i<nbErrors; ++i)
        {
            err=(*(modelFitter->GetError()))(i,0);
            if(fabs(err) > MAX_CALIB_ERR[calibIdx])
            {
                ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + CALIB_ERR_MSGE[calibIdx]);
            }
        }

        calib=calib->GetNextMethod();
        ++calibIdx;
    }
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CCSCalculator
///	Routine: GetOSWPortfolio
///	Returns: ARM_StdPortfolioPtr
///	Action : Get a diagonal swaption calibration portfolio
/////////////////////////////////////////////////////////////////
const ARM_StdPortfolioPtr ARM_CCSCalculator::GetOSWPortfolio(mdmKeysAlias oswModelKey) const
{
    /// Recall that Dom -> For -> Fx with "next" links
    switch(oswModelKey)
    {
    case OswDomModelKey:
        if(GetCalibMethod() != ARM_CalibMethodPtr(NULL))
            return GetCalibMethod()->GetPortfolio();

    case OswForModelKey:
        if( GetCalibMethod() != ARM_CalibMethodPtr(NULL) &&
            GetCalibMethod()->GetNextMethod() != NULL )
            return GetCalibMethod()->GetNextMethod()->GetPortfolio();
    }

    return ARM_StdPortfolioPtr(NULL);
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CCSCalculator
///	Routine: GetFxPortfolio
///	Returns: ARM_StdPortfolioPtr
///	Action : Get a fx option calibration portfolio
/////////////////////////////////////////////////////////////////
const ARM_StdPortfolioPtr ARM_CCSCalculator::GetFxPortfolio() const
{
    /// Recall that Dom -> For -> Fx with "next" links
    if( GetCalibMethod() != ARM_CalibMethodPtr(NULL) &&
        GetCalibMethod()->GetNextMethod() != NULL && 
        GetCalibMethod()->GetNextMethod()->GetNextMethod() != NULL ) 
        return GetCalibMethod()->GetNextMethod()->GetNextMethod()->GetPortfolio();
    else
        return ARM_StdPortfolioPtr(NULL);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CCSCalculator
///	Routine: GetFxCalib
///	Returns: ARM_CalibMethod*
///	Action : Get the fx calibration method
/////////////////////////////////////////////////////////////////
ARM_CalibMethod* ARM_CCSCalculator::GetFxCalib() const
{
    /// Recall that Dom -> For -> Fx with "next" links
    if( GetCalibMethod() != ARM_CalibMethodPtr(NULL) &&
        GetCalibMethod()->GetNextMethod() != NULL ) 
        return GetCalibMethod()->GetNextMethod()->GetNextMethod();
    else
        return NULL;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CCSCalculator
///	Routine: SetFxCalib
///	Returns: void
///	Action : Set the fx calibration method
/////////////////////////////////////////////////////////////////
void ARM_CCSCalculator::SetFxCalib(ARM_CalibMethod* fxCalib) const
{
    /// Recall that Dom -> For -> Fx with "next" links
    if( GetCalibMethod() != ARM_CalibMethodPtr(NULL) &&
        GetCalibMethod()->GetNextMethod() != NULL ) 
        GetCalibMethod()->GetNextMethod()->SetNextMethod(fxCalib);
    else
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : can't set the Fx calib method");
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

