/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file gencsocalculator.cpp
 *
 *  \brief file for the virtual calculator for Callable SpreadOption
 *
 *	\author  JP Riaudel
 *	\version 1.0
 *	\date July 2005
 */

#include "gpbase/removeidentifiedwarning.h"
#include "gpcalculators/gencsocalculator.h"



/// gpbase
#include "gpbase/env.h"
#include "gpbase/singleton.h"
#include "gpbase/autocleaner.h"
#include "gpbase/datestripcombiner.h"
#include "gpbase/gpmatrixlinalg.h"
#include "gpbase/globalconstant.h"
#include "gpbase/curveconvert.h"

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

/// kernel
#include "ccy/currency.h"
#include "inst/fixleg.h"
#include "inst/swapleg.h"
#include "mod/y2cmodel.h"
#include "inst/spreadoption.h"
//#include <inst/forex.h>
#include <util/fromto.h>

CC_BEGIN_NAMESPACE( ARM )

/// Reference schedules for CSO date structure
const unsigned int EXERCISE_SCHED =0;
const unsigned int NB_CSO_SCHED =1;

/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCSOCalculator
///	Routine: constructor
///	Returns: nothing
///	Action : build constructor
/////////////////////////////////////////////////////////////////
ARM_GenCSOCalculator::ARM_GenCSOCalculator(const ARM_Date& startDate,
										   const ARM_Date& endDate,
										   int CMSLong,
										   int CMSShort,
										   int cpnDayCount,
										   int cpnFreq,
										   int cpnResetTiming,
										   const ARM_Curve& cpnnominal,
										   const ARM_Curve& fixCoupon,
										   const ARM_Curve& leverageLong,
										   const ARM_Curve& leverageShort,
										   const ARM_Curve& cpnMin,
										   const ARM_Curve& cpnMax,
										   const ARM_Curve& strike,
										   int fundFreq,
										   int fundDaycount,
										   const ARM_Curve& fundSpread,
										   const ARM_Curve& fundLeverage,
										   int exerciseFreq,
										   int noticeGap,
										   int payRec,
										   const ARM_Curve& fees,
										   std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice,
										   vector<double> ModelDatas,
										   const ARM_MarketData_ManagerRep& mktDataManager,
										   const ARM_StringVector& mdmKeys)
: ARM_GenCalculator(mktDataManager),
	itsStartDate(startDate),
	itsFixEndDate(startDate),
	itsEndDate(endDate),	
	itsCMSLong(CMSLong),
	itsCMSShort(CMSShort),
	itsCpnDaycount(cpnDayCount),
	itsCpnFreq(cpnFreq),
	itsCpnResetTiming(cpnResetTiming),
	itsStubRule(K_SHORTSTART),
    itsCpnResetCal(""),
    itsCpnPayCal(""),
	itsCpnNominal(cpnnominal),
    itsFundNominal(cpnnominal),
	itsFixCpn(fixCoupon),
	itsLeverageLong(leverageLong),
	itsLeverageShort(leverageShort),
	itsCpnFloor(cpnMin),
	itsCpnCap(cpnMax),
	itsStrike(strike),
	itsFundFreq(fundFreq),
	itsFundDaycount(fundDaycount),
	itsFundSpread(fundSpread),
	itsFundLeverage(fundLeverage),
	itsExerciseFreq(exerciseFreq),
	itsNoticeGap(noticeGap),
	itsPayRec (payRec),
	itsFees(fees),
	itsProductsToPrice(productsToPrice),
	itsModelDatas(ModelDatas),
	itsHasBeenPriced(false),
	itsNbFixFlows(0),
	itsNbNoCall(0),
	itsFundingLeg(NULL),
	itsFixedLeg(NULL),
	itsInitFixedLeg(NULL),
	itsFloor(NULL),
	itsCap(NULL)
{
	SetName(ARM_CALLABLE_SPREADOPTION);

    /// Set keys for MDM datas access
	///  We suppose that OSW smiled model  and OSW Model are same
    if(mdmKeys.size() < NbKeys)
    {
        /// To be compatible with non basis version
        ARM_StringVector newMdMKeys(mdmKeys);
        newMdMKeys.resize(NbKeys);
        newMdMKeys[FundingKey]		  = newMdMKeys[YcKey];
        newMdMKeys[BasisKey]		  = newMdMKeys[YcKey];
		newMdMKeys[fundBasisKey]	  = newMdMKeys[YcKey];
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
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCSOCalculator
///	Routine: constructor (summit constructor in Fact used with InitFrom...)
///	Returns: nothing
///	Action : build calculator
///          without market data and model parameters
/////////////////////////////////////////////////////////////////
ARM_GenCSOCalculator::ARM_GenCSOCalculator(const ARM_Date& asOfDate,
										   const ARM_Date& startDate,
                                           const ARM_Date& fixEndDate,
										   const ARM_Date& endDate,
										   int CMSLong,
										   int CMSShort,
										   int cpnDayCount,
										   int cpnFreq,
										   int cpnResetTiming,
										   const ARM_Curve& cpnnominal,
										   const ARM_Curve& fixCoupon,
										   const ARM_Curve& leverageLong,
										   const ARM_Curve& leverageShort,
										   const ARM_Curve& cpnMin,
										   const ARM_Curve& cpnMax,
										   const ARM_Curve& strike,
										   int fundFreq,
										   int fundDaycount,
										   const ARM_Curve& fundSpread,
										   const ARM_Curve& fundLeverage,
										   int exerciseFreq,
										   int noticeGap,
										   int payRec,
										   const ARM_Curve& fees)
: ARM_GenCalculator(asOfDate),
	                        itsStartDate(startDate),
	                        itsFixEndDate(fixEndDate),
	                        itsEndDate(endDate),
	                        itsCMSLong(CMSLong),
	                        itsCMSShort(CMSShort),
	                        itsCpnDaycount(cpnDayCount),
	                        itsCpnFreq(cpnFreq),
	                        itsCpnResetTiming(cpnResetTiming),
	                        itsStubRule(K_SHORTSTART),
                            itsCpnResetCal(""),
                            itsCpnPayCal(""),
	                        itsCpnNominal(cpnnominal),
                            itsFundNominal(cpnnominal),
	                        itsFixCpn(fixCoupon),
	                        itsLeverageLong(leverageLong),
	                        itsLeverageShort(leverageShort),
	                        itsCpnFloor(cpnMin),
	                        itsCpnCap(cpnMax),
	                        itsStrike(strike),
	                        itsFundFreq(fundFreq),
	                        itsFundDaycount(fundDaycount),
	                        itsFundSpread(fundSpread),
	                        itsFundLeverage(fundLeverage),
	                        itsExerciseFreq(exerciseFreq),
	                        itsNoticeGap(noticeGap),
	                        itsPayRec (payRec),
	                        itsFees(fees),
	                        itsProductsToPrice(NULL),
	                        itsModelDatas(NULL),
	                        itsHasBeenPriced(false),
	                        itsNbFixFlows(0),
	                        itsNbNoCall(0),
	                        itsFundingLeg(NULL),
	                        itsFixedLeg(NULL),
	                        itsInitFixedLeg(NULL),
	                        itsFloor(NULL),
	                        itsCap(NULL)
{
	SetName(ARM_CALLABLE_SPREADOPTION);
/* TMP
    /// Schedules initialization
	DatesStructure();

    // fix period managing by fixEndDate
	while ( itsNbFixFlows < itsStructDateStrip->GetFlowEndDates()->size() 
            && 
			(*itsStructDateStrip->GetFlowEndDates())[itsNbFixFlows] < fixEndDate.GetJulian()+ ARM_GlobalConstant::ARM_SEVENDAYS_LAG )
		  itsNbFixFlows++;
*/
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCSOCalculator
///	Routine: constructor
///	Returns: nothing
///	Action : build constructor for basis case
/////////////////////////////////////////////////////////////////
ARM_GenCSOCalculator::ARM_GenCSOCalculator(const ARM_Date& startDate,
		const ARM_Date& fixEndDate,
		const ARM_Date& endDate,
		const ARM_Currency& CpnCcy,
		const ARM_Currency& FundCcy,
		int CMSLong,
		int CMSShort,
		int cpnDayCount,
		int cpnFreq,
		int cpnResetTiming,
		int stubRule,
		const string& cpnResetCal,
		const string& cpnPayCal,
		const ARM_Curve& cpnnominal,
		const ARM_Curve& leverageLong,
		const ARM_Curve& leverageShort,
		const ARM_Curve& cpnMin,
		const ARM_Curve& cpnMax,
		const ARM_Curve& strikes,
		int fundFreq,
		int fundDayCount,
		const ARM_Curve& fundnominal,
		const ARM_Curve& fundSpread,
		const ARM_Curve& fundLeverage,
		int exerciseFreq,
		int noticeGap,
		int payRec,
		size_t nbNCall,
		const ARM_Curve& fees,
		std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice,
		vector<double> ModelDatas,
		const ARM_MarketData_ManagerRep& mktDataManager)
:
    ARM_GenCalculator(mktDataManager.GetAsOfDate()),
		itsStartDate(startDate),
		itsEndDate(endDate),
		itsFixEndDate(fixEndDate),
		itsCMSLong(CMSLong),
		itsCMSShort(CMSShort),
		itsCpnDaycount(cpnDayCount),
		itsCpnFreq(cpnFreq),
		itsCpnResetTiming(cpnResetTiming),
		itsStubRule(stubRule),
		itsCpnResetCal(cpnResetCal),
		itsCpnPayCal(cpnPayCal),
		itsCpnNominal(cpnnominal),
		itsLeverageLong(leverageLong),
		itsLeverageShort(leverageShort),
		itsCpnFloor(cpnMin),
		itsCpnCap(cpnMax),
		itsStrike(strikes),
		itsFixCpn(strikes),
		itsFundFreq(fundFreq),
		itsFundDaycount(fundDayCount),
		itsFundNominal(fundnominal),
		itsFundSpread(fundSpread),
		itsFundLeverage(fundLeverage),
		itsExerciseFreq(exerciseFreq),
		itsNoticeGap(noticeGap),
		itsPayRec (payRec),
		itsFees(fees),
		itsProductsToPrice(productsToPrice),
		itsModelDatas(ModelDatas),
		itsHasBeenPriced(false),
		itsNbFixFlows(0),
		itsNbNoCall(nbNCall),
		itsFundingLeg(NULL),
		itsFixedLeg(NULL),
		itsInitFixedLeg(NULL),
		itsFloor(NULL),
		itsCap(NULL)
{
		SetName(ARM_CALLABLE_SPREADOPTION);

		/// nb fix flows
		for(int i(0); i <nbNCall; ++i) itsFees.GetOrdinates()[i] = NON_CALL_FEE;

		/// initialize and set MKt Keys
		ARM_StringVector keys(NbKeys,UNKNOWN_KEY_NAME);
		string cpnCcyName(CpnCcy.GetCcyName());
		string fundingCcyName(FundCcy.GetCcyName());
		keys[YcKey]			= YC_KEY_NAME + cpnCcyName;
		keys[MrsKey]		= MRS_KEY_NAME + cpnCcyName;
		keys[CorrelKey]		= CORREL_KEY_NAME + cpnCcyName;
        keys[CfModelKey]	= CFMODEL_KEY_NAME + cpnCcyName;
        keys[OswModelKey]	= OSWMODEL_KEY_NAME + cpnCcyName;
		keys[SoModelKey]	= SOMODEL_KEY_NAME + cpnCcyName;
		keys[VolRatioKey]	= VOLRATIO_KEY_NAME + cpnCcyName;
		keys[MrsSpreadKey]	= MRSSPREAD_KEY_NAME + cpnCcyName;
		
		///Set currencies
		SetCurrencyUnit(const_cast< ARM_Currency* >(&CpnCcy));
		SetFundingCcy(const_cast< ARM_Currency& >(FundCcy));

		if( fundingCcyName == cpnCcyName )
		{
			// no basis, no forex but keep compatibility !
			keys[FundingKey]	= keys[YcKey];
			keys[BasisKey]		= keys[YcKey];
			keys[fundBasisKey]  = keys[YcKey];
			// no basis, no forex
			SetForeignCcy(const_cast< ARM_Currency& >(CpnCcy));
			SetDomesticCcy(const_cast< ARM_Currency& >(CpnCcy));
		}
		else
		{
			keys[FundingKey]= YC_KEY_NAME + fundingCcyName;
			keys[BasisKey]= YC_BASIS_KEY_NAME + cpnCcyName;
			keys[fundBasisKey]= YC_BASIS_KEY_NAME + fundingCcyName;
			keys[ForexKey]= FOREX_KEY_NAME + cpnCcyName +"/"+ fundingCcyName;
			/// Set Foreign, basi and domestic Currencies
			ARM_Forex* forex = static_cast< ARM_Forex* >(mktDataManager.GetData(keys[ForexKey]));
			SetForeignCcy(*forex->GetMainCurrency());
			SetDomesticCcy(*forex->GetMoneyCurrency());
			ARM_Currency BasisCcy     = *static_cast< ARM_ZeroCurve* >(mktDataManager.GetData(keys[BasisKey]))->GetCurrencyUnit();
			SetBasisCcy(BasisCcy);
		} 
		SetKeys(keys);

		/// Schedules initialization
		DatesStructure();

		/// fix period managing by fixEndDate
		while(itsNbFixFlows < itsStructDateStrip->GetFlowEndDates()->size() && 
							   (*itsStructDateStrip->GetFlowEndDates())[itsNbFixFlows] < fixEndDate.GetJulian()+ ARM_GlobalConstant::ARM_SEVENDAYS_LAG )
			itsNbFixFlows++;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCSOCalculator
///	Routine: constructor
///	Returns: nothing
///	Action : build constructor for basis case
/////////////////////////////////////////////////////////////////
ARM_GenCSOCalculator::ARM_GenCSOCalculator(const ARM_Date& startDate,
		const ARM_Date& fixEndDate,
		const ARM_Date& endDate,
		const ARM_Currency& CpnCcy,
		const ARM_Currency& FundCcy,
		int CMSLong,
		int CMSShort,
		int cpnDayCount,
		int cpnFreq,
		int cpnResetTiming,
		int stubRule,
		const string& cpnResetCal,
		const string& cpnPayCal,
		const ARM_Curve& cpnnominal,
		const ARM_Curve& leverageLong,
		const ARM_Curve& leverageShort,
		const ARM_Curve& cpnMin,
		const ARM_Curve& cpnMax,
		const ARM_Curve& strikes,
		int fundFreq,
		int fundDayCount,
		const ARM_Curve& fundnominal,
		const ARM_Curve& fundSpread,
		const ARM_Curve& fundLeverage,
		int exerciseFreq,
		int noticeGap,
		int payRec,
		size_t nbNCall,
		const ARM_Curve& fees)
:
    ARM_GenCalculator(),
		itsStartDate(startDate),
		itsFixEndDate(startDate),
		itsEndDate(endDate),
		itsCMSLong(CMSLong),
		itsCMSShort(CMSShort),
		itsCpnDaycount(cpnDayCount),
		itsCpnFreq(cpnFreq),
		itsCpnResetTiming(cpnResetTiming),
		itsStubRule(stubRule),
		itsCpnResetCal(cpnResetCal),
		itsCpnPayCal(cpnPayCal),
		itsCpnNominal(cpnnominal),
		itsLeverageLong(leverageLong),
		itsLeverageShort(leverageShort),
		itsCpnFloor(cpnMin),
		itsCpnCap(cpnMax),
		itsStrike(strikes),
		itsFixCpn(*(ARM_Curve*)const_cast<ARM_Curve&>(strikes).Clone()),
		itsFundFreq(fundFreq),
		itsFundDaycount(fundDayCount),
		itsFundNominal(fundnominal),
		itsFundSpread(fundSpread),
		itsFundLeverage(fundLeverage),
		itsExerciseFreq(exerciseFreq),
		itsNoticeGap(noticeGap),
		itsPayRec (payRec),
		itsFees(fees),
		itsProductsToPrice(NULL),
		itsModelDatas(NULL),
		itsHasBeenPriced(false),
		itsNbFixFlows(0),
		itsFundingLeg(NULL),
		itsFixedLeg(NULL),
		itsInitFixedLeg(NULL),
		itsFloor(NULL),
		itsCap(NULL)
{
		SetName(ARM_CALLABLE_SPREADOPTION);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCSOCalculator
///	Routine: copy constructor
///	Returns: nothing
///	Action : copy constructor
/////////////////////////////////////////////////////////////////
ARM_GenCSOCalculator::ARM_GenCSOCalculator(const ARM_GenCSOCalculator& rhs)
: ARM_GenCalculator(rhs),
	itsStartDate(rhs.itsStartDate),
	itsFixEndDate(rhs.itsFixEndDate),
	itsEndDate(rhs.itsEndDate),
	itsCMSLong(rhs.itsCMSLong),
	itsCMSShort(rhs.itsCMSShort),
	itsCpnNominal(rhs.itsCpnNominal),
    itsFundNominal(rhs.itsFundNominal),
	itsFixCpn(rhs.itsFixCpn),  
	itsCpnDaycount(rhs.itsCpnDaycount),
	itsCpnFreq(rhs.itsCpnFreq),
	itsCpnResetTiming(rhs.itsCpnResetTiming),
	itsStubRule(rhs.itsStubRule),
	itsCpnFloor(rhs.itsCpnFloor),
	itsCpnCap(rhs.itsCpnCap),
	itsLeverageLong(rhs.itsLeverageLong),
	itsLeverageShort(rhs.itsLeverageShort),
	itsStrike(rhs.itsStrike),
	itsFundFreq(rhs.itsFundFreq),
	itsFundDaycount(rhs.itsFundDaycount),
	itsFundSpread(rhs.itsFundSpread),         
	itsFundLeverage(rhs.itsFundLeverage),         
    itsExerciseFreq(rhs.itsExerciseFreq),           
	itsNoticeGap(rhs.itsNoticeGap),
	itsFees(rhs.itsFees),
	itsFundDateStrip(rhs.itsFundDateStrip),         
    itsStructDateStrip(rhs.itsStructDateStrip),
	itsProductsToPrice(rhs.itsProductsToPrice),
	itsHasBeenPriced(rhs.itsHasBeenPriced),
	itsFirstEventIdx(rhs.itsFirstEventIdx),
	itsCSOPrice(rhs.itsCSOPrice),
	itsStructPrice(rhs.itsStructPrice),
	itsNbFixFlows(rhs.itsNbFixFlows),
	itsCalibDateStrip(rhs.itsCalibDateStrip),
	itsPayRec(rhs.itsPayRec),
	itsFundingLeg(CreateClone(rhs.itsFundingLeg)),
	itsFixedLeg(CreateClone(rhs.itsFixedLeg)),
	itsInitFixedLeg(CreateClone(rhs.itsInitFixedLeg)),
	itsFloor(CreateClone(rhs.itsFloor)),
	itsCap(CreateClone(rhs.itsCap))
{
}

////////////////////////////////////////////////////
///	Class  : ARM_GenCSOCalculator
///	Routine: operator =
///	Returns: 
///	Action : Assignment operator
////////////////////////////////////////////////////
ARM_GenCSOCalculator& ARM_GenCSOCalculator::operator = (const ARM_GenCSOCalculator& rhs)
{	
	if (&rhs != this)
	{ 
		this->~ARM_GenCSOCalculator();
        itsStartDate        = rhs.itsStartDate;
		itsFixEndDate       = rhs.itsFixEndDate;
	    itsEndDate          = rhs.itsEndDate;
	    itsCMSLong          = rhs.itsCMSLong;
	    itsCMSShort         = rhs.itsCMSShort;
	    itsCpnNominal       = rhs.itsCpnNominal;
        itsFundNominal      = rhs.itsFundNominal;
	    itsFixCpn           = rhs.itsFixCpn;  
	    itsCpnDaycount      = rhs.itsCpnDaycount;
	    itsCpnFreq          = rhs.itsCpnFreq;
		itsStubRule         = rhs.itsStubRule;
	    itsCpnFloor         = rhs.itsCpnFloor;
	    itsCpnCap           = rhs.itsCpnCap;
	    itsLeverageLong     = rhs.itsLeverageLong;
	    itsLeverageShort    = rhs.itsLeverageShort;
	    itsStrike           = rhs.itsStrike;
	    itsFundFreq         = rhs.itsFundFreq;
	    itsFundDaycount     = rhs.itsFundDaycount;
	    itsFundSpread       = rhs.itsFundSpread;         
		itsFundLeverage     = rhs.itsFundLeverage;         
        itsExerciseFreq     = rhs.itsExerciseFreq;           
	    itsNoticeGap        = rhs.itsNoticeGap;
	    itsFees             = rhs.itsFees;
	    itsFundDateStrip    = rhs.itsFundDateStrip;         
        itsStructDateStrip  = rhs.itsStructDateStrip;
	    itsProductsToPrice  = rhs.itsProductsToPrice;
	    itsHasBeenPriced    = rhs.itsHasBeenPriced;
	    itsFirstEventIdx    = rhs.itsFirstEventIdx;
	    itsCSOPrice         = rhs.itsCSOPrice;
	    itsStructPrice      = rhs.itsStructPrice;
	    itsNbFixFlows       = rhs.itsNbFixFlows;
	    itsCalibDateStrip   = rhs.itsCalibDateStrip;

        itsCpnModelKey      = rhs.itsCpnModelKey;
        itsFundingModelKey  = rhs.itsFundingModelKey;
        itsBasisRefModelKey = rhs.itsBasisRefModelKey;
		itsFixedLeg			= CreateClone(rhs.itsFixedLeg);
		itsFundingLeg		= CreateClone(rhs.itsFundingLeg);
		itsInitFixedLeg		= CreateClone(rhs.itsInitFixedLeg);
		itsFloor			= CreateClone(rhs.itsFloor);
		itsCap				= CreateClone(rhs.itsCap);
	}
	return *this;
}

ARM_GenCSOCalculator::~ARM_GenCSOCalculator() 
{
	delete  itsFundingLeg;
	itsFundingLeg = NULL;

	delete  itsFixedLeg;
	itsFixedLeg = NULL;

	delete  itsInitFixedLeg;
	itsInitFixedLeg = NULL;

	delete   itsFloor;
	itsFloor = NULL;

	delete   itsCap;
	itsCap = NULL;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCSOCalculator
///	Routine: SetModelKeys
///	Returns: void
///	Action : set the coupon, the funding & the basis reference
///          model names alias
/////////////////////////////////////////////////////////////////
void ARM_GenCSOCalculator::SetModelKeys()
{
    // Get the coupon & funding curve ccy
    string cpnCcyName(GetCurrencyUnit()->GetCcyName());
    string fundingCcyName(GetFundingCcy().GetCcyName());
    string basisCcyName(GetBasisCcy().GetCcyName());

    if(cpnCcyName == fundingCcyName)
    {
        // no basis effect
        itsBasisRefModelKey = YcKey;
        itsCpnModelKey      = YcKey;
        itsFundingModelKey  = YcKey;
    }
    else
    {
        if(cpnCcyName == basisCcyName)
        {
            itsBasisRefModelKey     = YcKey;        // Cpn forecast = YcKey (diffused)
            itsCpnModelKey          = BasisKey;     // Cpn discount = BasisKey
            itsFundingModelKey      = FundingKey;   // Funding forecast & discount = FundingKey
        }
        else
        {
            itsCpnModelKey          = YcKey;        // Cpn forecast & discount = YcKey (diffused)
            itsBasisRefModelKey     = FundingKey;   // Funding forecast = FundingKey
            itsFundingModelKey      = BasisKey;     // Funding discount = BasisKey
        }
    }
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCSOCalculator
///	Routine: CustomDatesStructure
///	Returns: ARM_DateStripCombiner
///	Action : create the list of all event dates of the CSO.
///			customized Schedule Not implemented yet.
/////////////////////////////////////////////////////////////////
ARM_DateStripCombiner ARM_GenCSOCalculator::CustomDatesStructure(const ARM_DateStripVector& dateStrips) const
{
	return ARM_DateStripCombiner();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCSOCalculator
///	Routine: DatesStructure
///	Returns: ARM_DateStripCombiner
///	Action : create the list of all event dates of the CSO.
/////////////////////////////////////////////////////////////////
ARM_DateStripCombiner ARM_GenCSOCalculator::DatesStructure() const
{
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	///
	/// generate schedules only once ...
	///
	if (itsExerciseDateStrip == ARM_DateStripPtr(NULL))
	{
		int fwdRule			= K_MOD_FOLLOWING;	// for forward dates
		int intRule			= K_ADJUSTED;			// for interest dates
		int resetTiming		= K_ADVANCE;
		int payTiming		= K_ARREARS;
		
		ARM_INDEX_TYPE indexType = ((ARM_Currency*)GetCurrencyUnit())->GetVanillaIndexType();
		char* DefaultresetCalendar = GetCurrencyUnit()->GetResetCalName(indexType);
		CC_NS(std,auto_ptr)<char> holdresetCalendar(DefaultresetCalendar);
		char* DefaultpayCalendar  = GetCurrencyUnit()->GetPayCalName(indexType);
		CC_NS(std,auto_ptr)<char> holdpayCalendar(DefaultpayCalendar);


		const char* resetCalendar = itsCpnResetCal == "" ? DefaultresetCalendar : itsCpnResetCal.c_str();
		const char* payCalendar	= itsCpnPayCal == "" ? DefaultpayCalendar: itsCpnPayCal.c_str();
		int resetGap			 = - GetCurrencyUnit()->GetSpotDays();

		/// exercise schedule
		ARM_DateStrip ExerSched(itsStartDate,itsEndDate,itsExerciseFreq,itsCpnDaycount,resetCalendar,fwdRule,intRule,
							    itsStubRule,itsNoticeGap,itsExerciseFreq, GETDEFAULTVALUE, payCalendar, resetTiming, payTiming);		
		
		/// past  no call
		size_t nbPastNoCall=0;
		std::vector<double>& exerciseDates = ExerSched.GetResetDates() ;
		while(nbPastNoCall < ExerSched.size() && (*exerciseDates)[nbPastNoCall] < asOfDate + K_NEW_DOUBLE_TOL)
			++nbPastNoCall;
		const_cast< ARM_GenCSOCalculator* >(this)->itsNbNoCall= itsNbNoCall > nbPastNoCall ? itsNbNoCall - nbPastNoCall : 0;

		/// Keep only the futur notification dates
		ExerSched.ResizeAndBuilt(nbPastNoCall,ExerSched.size());

		const_cast< ARM_GenCSOCalculator* >(this)->itsExerciseDateStrip = ARM_DateStripPtr(new ARM_DateStrip(ExerSched));		
		

		/// exercise schedule Unadj
		ARM_DateStrip ExerSchedUnadj(itsStartDate,itsEndDate,itsExerciseFreq,itsCpnDaycount,resetCalendar,fwdRule,K_UNADJUSTED,
							    itsStubRule,GETDEFAULTVALUE,itsExerciseFreq, GETDEFAULTVALUE, payCalendar, resetTiming, payTiming);		
		
		/// Keep only the futur notification dates
		ExerSchedUnadj.ResizeAndBuilt(nbPastNoCall,ExerSchedUnadj.size());

		const_cast< ARM_GenCSOCalculator* >(this)->itsExerciseDateStripUnadj = ARM_DateStripPtr(new ARM_DateStrip(ExerSchedUnadj));		
		

		/// funding schedule
		ARM_DateStrip FundingSched(	itsStartDate,itsEndDate,itsFundFreq,itsFundDaycount,resetCalendar,fwdRule,intRule,
								    itsStubRule,resetGap,itsFundFreq, GETDEFAULTVALUE, payCalendar, resetTiming, payTiming);

		nbPastNoCall=0;
		exerciseDates = FundingSched.GetResetDates() ;
		while(nbPastNoCall < ExerSched.size() && (*exerciseDates)[nbPastNoCall] < asOfDate + K_NEW_DOUBLE_TOL)
			++nbPastNoCall;
		FundingSched.ResizeAndBuilt(nbPastNoCall,FundingSched.size());

		const_cast< ARM_GenCSOCalculator* >(this)->itsFundDateStrip = ARM_DateStripPtr(new ARM_DateStrip(FundingSched));

		/// structured schedule
		ARM_DateStrip StructureSched(itsStartDate,itsEndDate,itsCpnFreq,itsCpnDaycount,resetCalendar,fwdRule,intRule,
									 itsStubRule,resetGap,itsCpnFreq, GETDEFAULTVALUE, payCalendar, resetTiming, payTiming);

		nbPastNoCall=0;
		exerciseDates = StructureSched.GetResetDates() ;
		while(nbPastNoCall < ExerSched.size() && (*exerciseDates)[nbPastNoCall] < asOfDate + K_NEW_DOUBLE_TOL)
			++nbPastNoCall;
		StructureSched.ResizeAndBuilt(nbPastNoCall,StructureSched.size());
		const_cast< ARM_GenCSOCalculator* >(this)->itsStructDateStrip = ARM_DateStripPtr( new ARM_DateStrip(StructureSched));

		/// calib sched
		const_cast< ARM_GenCSOCalculator* >(this)->itsCalibDateStrip = ARM_DateStripPtr( new ARM_DateStrip(StructureSched));

	}

    /// Merge schedules on "ResetDate"
    ARM_DateStripVector SchedVect(1);
    SchedVect[0]    = &*itsExerciseDateStrip;

    ARM_DateStripCombiner EventSchedule(SchedVect,"ResetDate");

	return EventSchedule;
}
/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCSOCalculator
///	Routine: CreateUnderlying
///	Returns: 
///	Action : Create a Power Reverse Swap
/////////////////////////////////////////////////////////////////
void ARM_GenCSOCalculator::CreateUnderlying()
{
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

	/// Funding Leg
	int indexType = !strcmp("EUR",GetFundingCcy().GetCcyName()) ? FromFrequencyToEuriborType(itsFundFreq) :
																 FromFrequencyToLiborType(itsFundFreq);
	
	ARM_ReferenceValue  rspread(GPCurveToRefValue(itsFundSpread,asOfDate,100.0));
	ARM_ReferenceValue  fundNominal(GPCurveToRefValue(itsFundNominal,asOfDate));

	delete itsFundingLeg;
	itsFundingLeg = new ARM_SwapLeg(itsStartDate,
							itsEndDate,
							(ARM_INDEX_TYPE)indexType,
							K_RCV,
							0.0,
							itsFundFreq,
							itsFundFreq,
							K_ADVANCE,
							K_ARREARS,
							&GetFundingCcy(),
							K_ADJUSTED,
							10000,
							const_cast<char*>(itsCpnResetCal.c_str()),
							const_cast<char*>(itsCpnPayCal.c_str()),
							1,
							K_NX_NONE,
							itsStubRule);	
	itsFundingLeg->SetVariableSpread(&rspread);	
	itsFundingLeg->SetAmount(&fundNominal);

	ARM_ReferenceValue  fixedRates(GPCurveToRefValue(itsStrike,asOfDate,100.0));
	ARM_ReferenceValue  nominal(GPCurveToRefValue(itsCpnNominal,asOfDate));

	delete itsInitFixedLeg;
	itsInitFixedLeg =  new ARM_FixLeg(itsStartDate,
				itsFixEndDate, 
                &fixedRates,
				K_RCV, 
				itsCpnFreq,
				itsCpnDaycount,
				K_COMP_PROP,
				K_ARREARS,
				K_ADJUSTED,
				itsStubRule,
				&GetDomesticCcy(),
				const_cast<char*>(itsCpnPayCal.c_str()),
				K_NX_NONE,
				"NULL",
				1);
	((ARM_FixLeg*)itsInitFixedLeg)->SetVarCoupons(&fixedRates);
	itsInitFixedLeg->SetAmount(&nominal);

	ARM_ReferenceValue  refFloor(GPCurveToRefValue(itsCpnFloor,asOfDate,100.0));

	delete itsFixedLeg;
	itsFixedLeg =   new ARM_FixLeg(itsFixEndDate,
				itsEndDate, 
                &refFloor,
				K_RCV, 
				itsCpnFreq,
				itsCpnDaycount,
				K_COMP_PROP,
				K_ARREARS,
				K_ADJUSTED,
				itsStubRule,
				&GetDomesticCcy(),
				const_cast<char*>(itsCpnPayCal.c_str()),
				K_NX_NONE,
				"NULL",
				1);
	itsFixedLeg->SetAmount(&nominal);

	ARM_Vector* fixing1 = NULL;
	ARM_Vector* fixing2 = NULL;

	ARM_ReferenceValue  shortlvge(GPCurveToRefValue(itsLeverageShort,asOfDate));
	ARM_ReferenceValue  longlvge(GPCurveToRefValue(itsLeverageLong,asOfDate));

	ARM_Curve cpnFloor(itsCpnFloor);	
	cpnFloor.GetOrdinates()+= itsStrike.GetOrdinates();

	ARM_ReferenceValue  refcpnFloor(GPCurveToRefValue(cpnFloor,asOfDate,100.0));
	
	delete itsFloor;
	itsFloor = new ARM_SpreadOption( itsFixEndDate,
			   itsEndDate,
			   K_CAP,
			   &refcpnFloor,
			   (ARM_INDEX_TYPE)itsCMSShort,
			   (ARM_INDEX_TYPE)itsCMSLong,
			   &shortlvge,
			   &longlvge, 
			   itsCpnDaycount,
			   itsCpnFreq,
			   itsCpnFreq,
			   itsCpnResetTiming,
			   K_ARREARS,
			   GetCurrencyUnit(),
			   1.0,
			   fixing1,
			   fixing2);
	itsFloor->SetAmount(&nominal);

	ARM_Curve cpnCap(itsCpnCap);	
	cpnCap.GetOrdinates()+= itsStrike.GetOrdinates();

    ARM_ReferenceValue  refcpnCap(GPCurveToRefValue(cpnCap,asOfDate,100.0));

	delete itsCap;
	itsCap = new ARM_SpreadOption ( itsFixEndDate,
				itsEndDate,
				K_CAP,
			   &refcpnCap,
			   (ARM_INDEX_TYPE)itsCMSShort,
			   (ARM_INDEX_TYPE)itsCMSLong,
			   &shortlvge,
			   &longlvge, 
			   itsCpnDaycount,
			   itsCpnFreq,
			   itsCpnFreq,
			   itsCpnResetTiming,
			   K_ARREARS,
			   GetCurrencyUnit(),
			   1.0,
			   fixing1,
			   fixing2);
	itsCap->SetAmount(&nominal);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCSOCalculator
///	Routine: DatesStructure
///	Returns: ARM_DateStripCombiner
///	Action : create the list of all event dates of the CSO.
/////////////////////////////////////////////////////////////////
ARM_Vector* ARM_GenCSOCalculator::ComputeAll()
{
	ARM_Vector* result = new ARM_Vector(5);
	CreateUnderlying();
	ARM_ZeroCurve* fundzcCurve = static_cast< ARM_ZeroCurve* >(GetMktDataManager()->GetData(GetKeys()[FundingKey]));
	ARM_ZeroCurve* fundBsZcCurve = static_cast< ARM_ZeroCurve* >(GetMktDataManager()->GetData(GetKeys()[fundBasisKey]));

	ARM_Y2CModel ycModel(fundzcCurve,fundBsZcCurve);
	itsFundingLeg->SetModel(&ycModel);
	double fundPrice = itsFundingLeg->ComputePrice();
	if(IsBasis()){
		double spot = static_cast< ARM_Forex* >(GetMktDataManager()->GetData(GetKeys()[ForexKey]))->GetMarketPrice();
		double ZcStart = fundBsZcCurve->DiscountPrice(itsStartDate);
		double ZcEnd = fundBsZcCurve->DiscountPrice(itsEndDate);
		double notional = itsFundNominal.GetOrdinates()[0];
		double NxAtBegin = notional*ZcStart;
		double NxAtEnd = notional*ZcEnd;
		fundPrice += (NxAtEnd - NxAtBegin);
		fundPrice *= spot;

		ARM_ZeroCurve* BsZcCurve = static_cast< ARM_ZeroCurve* >(GetMktDataManager()->GetData(GetKeys()[BasisKey])); 
		ZcStart = BsZcCurve->DiscountPrice(itsStartDate);
		ZcEnd = BsZcCurve->DiscountPrice(itsEndDate);
		notional = itsCpnNominal.GetOrdinates()[0];
		NxAtBegin = notional*ZcStart;
		NxAtEnd = notional*ZcEnd;

		fundPrice += (NxAtBegin-NxAtEnd);
	}
	(*result)[1] = fundPrice;

	ARM_BSModel* SOBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[SoModelKey]) );

	itsInitFixedLeg->SetModel(SOBSModel);
	double initPrice = itsInitFixedLeg->ComputePrice();

	itsFixedLeg->SetModel(SOBSModel);
	double fixedPrice = itsFixedLeg->ComputePrice();
	(*result)[2] = fixedPrice + initPrice;

	itsFloor->SetModel(SOBSModel);
	double floorPrice = itsFloor->ComputePrice();
	(*result)[3] = floorPrice;

	itsCap->SetModel(SOBSModel);
	double capPrice = itsCap->ComputePrice();
	(*result)[4] = capPrice;

	double underlyingPrice = fundPrice - (fixedPrice + initPrice + floorPrice - capPrice);

	(*result)[0] = underlyingPrice;

	return result;
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
