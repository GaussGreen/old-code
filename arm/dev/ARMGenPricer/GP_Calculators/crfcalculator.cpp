/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file crfcalculator.cpp
 *
 *  \brief file for the CRF calculator
 *	\author  JM Prie
 *	\version 1.0
 *	\date March 2004
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpcalculators/crfcalculator.h"

/// gpcalculators
#include "gpcalculators/argconvdefault.h"

/// gpbase
#include "gpbase/autocleaner.h"
#include "gpbase/singleton.h"
#include "gpbase/ostringstream.h"
#include "gpbase/numericconstant.h"
#include "gpbase/curve.h"
#include "gpbase/datestripcombiner.h"
#include "gpbase/datestrip.h"
#include "gpbase/env.h"
#include "gpbase/gplinalgconvert.h"
#include "gpbase/surface.h"
#include "gpbase/surfacetypedef.h"
#include "gpbase/cloneutilityfunc.h"
#include "gpbase/utilityport.h"
#include "gpbase/globalconstant.h"


/// gpinfra
#include "gpinfra/mktdatamanagerrep.h"
#include "gpinfra/argconvdefault.h"
#include "gpinfra/gensecurity.h"
#include "gpinfra/dealdescription.h"
#include "gpinfra/pricingadviser.h"
#include "gpinfra/modelnrefcall.h"
#include "gpinfra/gramnode.h"
#include "gpinfra/genpricer.h"
#include "gpinfra/numeraire.h"
#include "gpinfra/numerairefactory.h"
#include "gpinfra/cstmanager.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/surfacemodelparam.h"
#include "gpinfra/modelnamemap.h"
#include "gpinfra/pricerinfo.h"


/// gpcalib
#include "gpcalib/calibmethod.h"
#include "gpcalib/vanillacap.h"
#include "gpcalib/kerneltogp.h"
#include "gpcalib/modelfitterdes.h"
#include "gpcalib/vanillapricer.h"
#include "gpcalib/bootstrap1d.h"

/// gpmodels
#include "gpmodels/multiassets.h"
#include "gpmodels/hw1f.h"
#include "gpmodels/hwfactory.h"
#include "gpmodels/modelparamshw1f.h"
#include "gpmodels/QGM1F.h"
#include "gpmodels/Normal_Model.h"
#include "gpmodels/Normal_ModelParams.h"
#include "gpmodels/ModelParamsQGM1F.h"
#include "gpmodels/ForwardMarginBasis.h"
#include "gpmodels/ForwardMarginIR.h"
#include "gpmodels/forwardforex.h"
#include "gpmodels/hybridbasisfwdir.h"
#include "gpmodels/sfrm.h"
#include "gpmodels/modelparamssfrm.h"
#include "gpmodels/modelparamssfrmfactory.h"
#include "gpmodels/MarketIRModel.h"

/// gpnummethods
#include "gpnummethods/treemethod.h"
#include "gpnummethods/treebase.h"
#include "gpnummethods/treefactory.h"
#include "gpnummethods/scheduler.h"
#include "gpnummethods/sampler.h"
#include "gpnummethods/truncator.h"
#include "gpnummethods/reconnector.h"
#include "gpnummethods/smoother.h"

/// gpnumlib
#include "gpnumlib/solver.h"
#include "gpnumlib/brent.h"

/// gpclosedforms
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/qgm_skewcalibration.h"

/// kernel
#include <inst/swaption.h>
#include <inst/fixleg.h>
#include <ccy/currency.h>
#include <crv/volflat.h>
#include <mod/bssmiled.h>
#include <inst/forex.h>
#include <glob/paramview.h>


/// STL
#include <iomanip> /// for setprecision()
#include <list>
CC_USING_NS(std,list)
#include <map>
CC_USING_NS( std, map )
#include <memory>


CC_BEGIN_NAMESPACE( ARM )


/// To enable/disable equivalent strike computation for
/// diagonal swaptions with moment matching for 2 LN
#define IS_LN_APPROX_K_EQUIV

/// To define a non callable event date
const double NON_CALL_FEE=1.0e+15;

/// To indicate that exercise gap is not input => callable
/// event dates must be inferred from a given exercise schedule
const int NOT_DEF_EXER_GAP=-1111;

/// To indicate no link to any additional standard swaps for
/// equivalent strike/vol computation (diagonal swaption generation)
const int NO_OTHER_STDSWAPS=-1;


/// Reference schedules for CRF date structure
const unsigned int FIX_CPN_SCHED=0;
const unsigned int RF_CPN_SCHED =1;
const unsigned int EXER_SCHED   =2;
const unsigned int NB_CRF_SCHED =3;

/// Type of CRF event
const unsigned int FIX_CPN_EVENT  = 1;
const unsigned int RF_CPN_EVENT   = 2;
const unsigned int EXER_EVENT     = 4;
const unsigned int CPN_EVENT      = FIX_CPN_EVENT | RF_CPN_EVENT;


/// Tree using by default 30 time steps per year
const int TREE_NBSTEPS_PER_YEAR=30;


/// H&W & QGM sigma range [10bp,500bp] with a 100bp default value
/// At initialisation step, vols will be set in range [50bp,150bp]
const double SIGMA_LOWER_BOUND      = 0.001;
const double SIGMA_UPPER_BOUND      = 0.05;
const double SIGMA_DEFAULT_VALUE    = 0.01;
const double SIGMA_MIN_INIT_VALUE   = 5.0e-4;
const double SIGMA_MAX_INIT_VALUE   = 250.0e-4;

/// H&W & QGM MRS range [-15%,50%] with a 0% default value
const double MRS_LOWER_BOUND        = -0.15;
const double MRS_UPPER_BOUND        = 0.5;
const double LEVEL_TO_CHOOSE_PORTFOLIO = 1.0;

const double MRS_DEFAULT_VALUE      = 0.0;

/// QGM Skew range [-15,40] with a 0 (<=>H&W) default value
const double SKEW_LOWER_BOUND      = 1.0e-10;
const double SKEW_UPPER_BOUND      = 150.0;
const double SKEW_DEFAULT_VALUE    = 25.0;

/// QGM Default Parameters for the calibration
const double SKEW_CALIB_PRECISION   = 1.0e-6;
const double VOL_CALIB_PRECISION    = 1.0e-6;
const size_t MAXITER				= 30;
const size_t NB_ITER_MAX			= 5;

/// Vega limit to be selected in portfolio for volatility bootstrapping 
const double VEGA_MIN_TO_SELECT=1.0e-5; 
const double VEGA_MIN_TO_REMOVE=1.0e-6; 
const double OSW_DEFAULT_WEIGHT=1.0;
const double OSW_DEFAULT_PRICE=1.0e+100;

/// 0.5bp (expected 0.01bp) of kappa to be selected in portfolio for strike spread adjustments
const double KAPPA_MIN_TO_SELECT=0.00005; //0.000001; => not activated due to regressions...
const double CF_DEFAULT_WEIGHT=1.0;


/// Type of product able to be priced by the CRF deal description
const unsigned int CRF_PRICE        = 0;
const unsigned int CAP_PRICE        = 1;
const unsigned int FLOOR_PRICE      = 2;
const unsigned int FUNDING_PRICE    = 3;
const unsigned int STDLEG_PRICE     = 4;
const unsigned int RFLEG_PRICE      = 5;
const unsigned int STDSWAP_PRICE    = 6;
const unsigned int RFSWAP_PRICE     = 7;
const unsigned int BERMUDA_PRICE    = 8;


/// Equivalent datas computed for diagonal swaption
const unsigned int OSW_SWAP_RATE    = 0;
const unsigned int OSW_TARGET_VOL   = 1;
const unsigned int OSW_TARGET_PRICE = 2;
const unsigned int OSW_TARGET_VEGA  = 3;
const unsigned int OSW_NB_EQUIVDATA = 4;


/// Default MDM key names
const string YC_KEY_NAME            = "YC_";
const string OSWMODEL_KEY_NAME      = "OSWMOD_";
const string CFMODEL_KEY_NAME       = "CFMOD_";
const string MRS_KEY_NAME           = "MRS_";
const string YC_BASIS_KEY_NAME      = "YC_BASIS_";
const string FOREX_KEY_NAME         = "FOREX_";
const string UNKNOWN_KEY_NAME       = "UNKNOWN";


/////////////////////////////////////////////////////////////////
///	Class  : CRFStrikeAdjFct
///	Routine: constructor & destructor
///	Returns: nothing
///	Action : 
/////////////////////////////////////////////////////////////////
CRFStrikeAdjFct::CRFStrikeAdjFct(size_t prodIdx, ARM_CRFCalculator* crfCalculator)
: itsProdIdx(prodIdx), itsCRFCalculator(crfCalculator)
{}

CRFStrikeAdjFct::~CRFStrikeAdjFct()
{}


/////////////////////////////////////////////////////////////////
///	Class  : CRFStrikeAdjFct
///	Routine: () operator
///	Returns: a double
///	Action : compute the error between the model & market prices
///          of a current caplet or floorlet
/////////////////////////////////////////////////////////////////
double CRFStrikeAdjFct::operator () ( double x ) const
{
    /// Set the strike
    ARM_VanillaCapArg* cf = static_cast< ARM_VanillaCapArg* >( &(*(itsCRFCalculator->GetVanillaArgVect()[itsProdIdx])) );
	(*(cf->GetStrikes()))[0]=x;

    double modelPrice = itsCRFCalculator->GetVanillaArgVect()[itsProdIdx]->Price(&(*(itsCRFCalculator->GetPricingModel())));
    double mktPrice = itsCRFCalculator->GetVanillaArgVect()[itsProdIdx]->GetMktPrice();

    return modelPrice - mktPrice;
}



const string ARM_CRFCalculator::CRFColNamesTable [] =
{
    "EventDate",
    "StartDate",
    "EndDate",
    "PayDate",
    "IT",
    "IndexStartDate",
    "Strike",
    "CpnMin",
    "CpnMax",
    "Leverage",
    "AdjCapSpread",
    "AdjFloorSpread",
    "Nominal",
    "FundingSpread",
    "FundingNominal",
    "FundingStartDate",
    "FundingEndDate",
    "FundingAnnuity",
    "FundingNominalExchange",
    "FundingVarFlow",
    "FundingSpreadFlow",
    "Funding",
    "DFPay",
    "CpnIndex",
    "CpnIndexFlow",
    "StrikeFlow",
    "StdVarFlow",
    "StdFixFlow",
    "Caplet",
    "Floorlet",
    "StdFlow",
    "RFFlow",
    "StdSwaplet",
    "StdSwap",
    "RFSwaplet",
    "RFSwap",
    "Fee",
    "StdBermuda",
	"FinalDate",
	"SwapRate",	
    "RFBermuda",
	"Frontier",
    "StdBermudaPrice",
    "RFBermudaPrice"
};


////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: Copy constuctor
///	Returns: void
///	Action : Arguments copy
////////////////////////////////////////////////////
ARM_CRFCalculator::ARM_CRFCalculator(const ARM_CRFCalculator& rhs)
:	ARM_GenCalculator( rhs ),
	itsStartDate( rhs.itsStartDate),
	itsEndDate( rhs.itsEndDate),
    itsNominal( rhs.itsNominal),
    itsFundNominal( rhs.itsFundNominal),
    itsStrike( rhs.itsStrike),
    itsPayRec( rhs.itsPayRec),
    itsFixEndDate( rhs.itsFixEndDate),
    itsFixDayCount( rhs.itsFixDayCount),
    itsCpnDayCount( rhs.itsCpnDayCount),
    itsCpnFreq( rhs.itsCpnFreq),
    itsCpnTiming( rhs.itsCpnTiming),
    itsCpnIndexTerm( rhs.itsCpnIndexTerm),
    itsCpnIndexDayCount( rhs.itsCpnIndexDayCount),
    itsCpnResetCal( rhs.itsCpnResetCal),
    itsCpnPayCal( rhs.itsCpnPayCal),
	itsStubRule	( rhs.itsStubRule),
    itsCpnResetGap( rhs.itsCpnResetGap),
    itsLeverage( rhs.itsLeverage),
    itsCpnMin( rhs.itsCpnMin),
    itsCpnMax( rhs.itsCpnMax),
    itsFundSpread( rhs.itsFundSpread),
    itsFundFreq( rhs.itsFundFreq),
    itsFundDayCount	( rhs.itsFundDayCount),
    itsExerGap( rhs.itsExerGap),
    itsNbNonCall( rhs.itsNbNonCall),
    itsExerFee( rhs.itsExerFee),
	itsExerCal( rhs.itsExerCal),
	itsModelType( rhs.itsModelType),
	itsMRSStrikeType( rhs.itsMRSStrikeType),
	itsSkewCalFlag( rhs.itsSkewCalFlag),
	itsSkewReCalFlag(rhs.itsSkewReCalFlag),
	itsIsFrontier(rhs.itsIsFrontier),
	 itsProductToPrice(rhs.itsProductToPrice),
	itsOSWCalibFlag( rhs.itsOSWCalibFlag),
    itsMRSCalibType	( rhs.itsMRSCalibType),
    itsCapCalibFlag( rhs.itsCapCalibFlag),
    itsFloorCalibFlag( rhs.itsFloorCalibFlag),
    itsCpnModelKey(rhs.itsCpnModelKey),
    itsFundingModelKey( rhs.itsFundingModelKey),
    itsBasisRefModelKey (rhs.itsBasisRefModelKey),
    itsFirstCallIdx(rhs.itsFirstCallIdx),
	 itsOtherStdSwapsIdx(rhs.itsOtherStdSwapsIdx),
	itsKeepOSWIdx(rhs.itsKeepOSWIdx),
	itsvStrike(rhs.itsvStrike),
	itsCapFloorPF(CreateClone(rhs.itsCapFloorPF)),
	itsOWSSkewPF(CreateClone(rhs.itsOWSSkewPF)) 
{

	itsConvexityModel = ARM_PricingModelPtr( rhs.itsConvexityModel != ARM_PricingModelPtr(NULL) ? static_cast< ARM_PricingModel* >(rhs.itsConvexityModel->Clone()) : NULL );

    itsVanillaArgVect.resize(rhs.itsVanillaArgVect.size());
    for(int i=0;i<rhs.itsVanillaArgVect.size();++i)
        itsVanillaArgVect[i] = static_cast< ARM_VanillaArg* >(const_cast< ARM_CRFCalculator& >(rhs).itsVanillaArgVect[i]->Clone());

    itsStdSwaps.resize(rhs.itsStdSwaps.size());
    for(i=0;i<rhs.itsStdSwaps.size();++i)
    {
        itsStdSwaps[i].first=static_cast< ARM_Swap* >(const_cast< ARM_CRFCalculator& >(rhs).itsStdSwaps[i].first->Clone());
        itsStdSwaps[i].second=rhs.itsStdSwaps[i].second;
    }

    itsFundStdSwaps.resize(rhs.itsFundStdSwaps.size());
    for(i=0;i<rhs.itsFundStdSwaps.size();++i)
    {
        itsFundStdSwaps[i].first=static_cast< ARM_Swap* >(const_cast< ARM_CRFCalculator& >(rhs).itsFundStdSwaps[i].first->Clone());
        itsFundStdSwaps[i].second=rhs.itsFundStdSwaps[i].second;
    }

    itsIndexStdSwaps.resize(rhs.itsIndexStdSwaps.size());
    for(i=0;i<rhs.itsIndexStdSwaps.size();++i)
    {
        itsIndexStdSwaps[i].first=static_cast< ARM_Swap* >(const_cast< ARM_CRFCalculator& >(rhs).itsIndexStdSwaps[i].first->Clone());
        itsIndexStdSwaps[i].second=rhs.itsIndexStdSwaps[i].second;
    }
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: Constructor
///	Returns: void
///	Action : builds the object (stand alone version)
/////////////////////////////////////////////////////////////////
ARM_CRFCalculator::ARM_CRFCalculator(const ARM_Date& startDate,
                const ARM_Date& endDate,
                const ARM_ReferenceValue& strike,
                int payRec,
                const ARM_Date& fixEndDate,
                int fixDayCount,
                int cpnDayCount,
                int cpnFreq,
                int cpnTiming,
                const string& cpnIndexTerm,
                int cpnIndexDayCount,
                const string& cpnResetCal,
                const string& cpnPayCal,
				int stubRule,
                int cpnResetGap,
                const ARM_ReferenceValue& leverage,
                const ARM_ReferenceValue& cpnMin,
                const ARM_ReferenceValue& cpnMax,
                const ARM_ReferenceValue& fundSpread,
                int fundFreq,
                int fundDayCount,
                const ARM_ReferenceValue& nominal,
                int exerGap,
                int nbNonCall,
                const ARM_ReferenceValue& exerFee,
                ARM_SigmaCalibType oswCalibType,
                ARM_MRSCalibType mrscalibType,
				ARM_MRSStrikeCalibType mrsStrikeType,
                bool capCalibFlag,
                bool floorCalibFlag,
                const ARM_MarketData_ManagerRep& mktDataManager,
                const ARM_StringVector& mdmKeys,
                const ARM_ReferenceValue& fundnominal,
				ARM_ModelType modelType,
				bool skewCalFlag,
				long isFrontier)
:	ARM_GenCalculator(mktDataManager),
	itsStartDate(startDate),
	itsEndDate(endDate),
    itsStrike(strike),
    itsPayRec(payRec),
    itsFixEndDate(fixEndDate),
    itsFixDayCount(fixDayCount),
    itsCpnDayCount(cpnDayCount),
    itsCpnFreq(cpnFreq),
    itsCpnTiming(cpnTiming),
    itsCpnIndexTerm(cpnIndexTerm),
    itsCpnIndexDayCount(cpnIndexDayCount),
    itsCpnResetCal(cpnResetCal),
    itsCpnPayCal(cpnPayCal),
	itsStubRule(stubRule),
    itsCpnResetGap(cpnResetGap),
    itsLeverage(leverage),
    itsCpnMin(cpnMin),
    itsCpnMax(cpnMax),
    itsFundSpread(fundSpread),
    itsFundFreq(fundFreq),
    itsFundDayCount(fundDayCount),
    itsNominal(nominal),
    itsFundNominal(fundnominal),
    itsExerGap(exerGap),
    itsNbNonCall(nbNonCall),
    itsExerFee(exerFee),
    itsOSWCalibFlag(oswCalibType),
    itsMRSCalibType(mrscalibType),
    itsCapCalibFlag(capCalibFlag),
    itsFloorCalibFlag(floorCalibFlag),
	itsModelType(modelType),
	itsMRSStrikeType(mrsStrikeType),
	itsSkewCalFlag(skewCalFlag),
	itsSkewReCalFlag(false),
	itsIsFrontier(isFrontier),
    itsProductToPrice(CRF_PRICE),
    itsCapFloorPF(NULL),
	itsOWSSkewPF(NULL),
    itsFirstCallIdx(0),
	itsExerCal(cpnResetCal),
	itsKeepOSWIdx(),
	itsvStrike(),
	itsVanillaArgVect(0)
{
    SetName(ARM_CALLREVFLOATER);

    /// Set keys for MDM datas access
    if(mdmKeys.size() < NbKeys)
    {
        /// To be compatible with non basis version
        ARM_StringVector newMdMKeys(mdmKeys);
        newMdMKeys.resize(NbKeys);
        newMdMKeys[FundingKey]  = newMdMKeys[YcKey];
        newMdMKeys[BasisKey]    = newMdMKeys[YcKey];
        newMdMKeys[ForexKey]    = UNKNOWN_KEY_NAME;
        SetKeys(newMdMKeys);
    }
    else
        SetKeys(mdmKeys);
    
    /// Check input datas
    CheckDataAndTimeIt();

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

    /// Set the coupon & funding model name alias i.e. keys to access to
    /// the right model for each GP keyword (1st argument)
    SetModelKeys();

    /// Create the Generic Security paid in coupon currency
    CreateAndSetDealDescriptionAndTimeIt(GetKeys()[itsCpnModelKey]);

	/// Create and initialize the pricing model
	CreateAndSetModelAndTimeIt();
	/// Create the calibration set for volatility bootstapping and
	/// strike spread adjustments
	CreateAndSetCalibrationAndTimeIt();
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: Constructor
///	Returns: void
///	Action : builds the object (generic security only)
/////////////////////////////////////////////////////////////////
ARM_CRFCalculator::ARM_CRFCalculator(const ARM_Date& asOfDate,
                                        const ARM_Date& startDate,
                                        const ARM_Date& endDate,
                                        const ARM_ReferenceValue& strike,
                                        int payRec,
                                        const ARM_Date& fixEndDate,
                                        int fixDayCount,
                                        int cpnDayCount,
                                        int cpnFreq,
                                        int cpnTiming,
                                        const string& cpnIndexTerm,
                                        int cpnIndexDayCount,
                                        const string& cpnResetCal,
                                        const string& cpnPayCal,
				                        int stubRule,
                                        int cpnResetGap,
                                        const ARM_ReferenceValue& leverage,
                                        const ARM_ReferenceValue& cpnMin,
                                        const ARM_ReferenceValue& cpnMax,
                                        const ARM_ReferenceValue& fundSpread,
                                        int fundFreq,
                                        int fundDayCount,
                                        const ARM_ReferenceValue& nominal,
                                        const ARM_ReferenceValue& exerFee,
                                        const ARM_Currency& cpnCcy,     
                                        const ARM_Currency& fundingCcy, 
                                        const ARM_Currency& basisCcy,   
                                        const ARM_Currency& domesticCcy,
                                        const ARM_Currency& foreignCcy, 
                                        const ARM_ReferenceValue& fundnominal,
				                        ARM_ModelType modelType,
                                        bool SkewReCalFlag)
                        :	ARM_GenCalculator(asOfDate),
	                        itsStartDate(startDate),
	                        itsEndDate(endDate),
                            itsStrike(strike),
                            itsPayRec(payRec),
                            itsFixEndDate(fixEndDate),
                            itsFixDayCount(fixDayCount),
                            itsCpnDayCount(cpnDayCount),
                            itsCpnFreq(cpnFreq),
                            itsCpnTiming(cpnTiming),
                            itsCpnIndexTerm(cpnIndexTerm),
                            itsCpnIndexDayCount(cpnIndexDayCount),
                            itsCpnResetCal(cpnResetCal),
                            itsCpnPayCal(cpnResetCal),
	                        itsStubRule(stubRule),
                            itsCpnResetGap(cpnResetGap),
                            itsLeverage(leverage),
                            itsCpnMin(cpnMin),
                            itsCpnMax(cpnMax),
                            itsFundSpread(fundSpread),
                            itsFundFreq(fundFreq),
                            itsFundDayCount(fundDayCount),
                            itsNominal(nominal),
                            itsFundNominal(fundnominal),
                            itsExerGap(NOT_DEF_EXER_GAP),
                            itsNbNonCall(NOT_DEF_EXER_GAP),
                            itsExerFee(exerFee),
                            itsOSWCalibFlag(ARM_SigmaCalibrationType::strikeEquivalent),
                            itsMRSCalibType(ARM_MRSCalibrationType::Unknown),
                            itsCapCalibFlag(true),
                            itsFloorCalibFlag(false),
                            itsProductToPrice(CRF_PRICE),
	                        itsCapFloorPF(NULL),
                            itsFirstCallIdx(0),
	                        itsOWSSkewPF(NULL),
                            itsModelType(modelType),
	                        itsMRSStrikeType(ARM_MRSStrikeCalibrationType::strikeEquivalent),
	                        itsSkewCalFlag(true),
	                        itsSkewReCalFlag(SkewReCalFlag),
							itsIsFrontier(0),
	                        itsExerCal(cpnResetCal),
	                        itsKeepOSWIdx(),
							itsvStrike()
{
    SetName(ARM_CALLREVFLOATER);

    /// Set the coupon/payment currency (inherited from ARM_Security)
    SetCurrencyUnit(const_cast< ARM_Currency* >(&cpnCcy));
    SetFundingCcy(const_cast< ARM_Currency& >(fundingCcy));
    SetBasisCcy(const_cast< ARM_Currency& >( basisCcy));
    SetDomesticCcy(const_cast< ARM_Currency& > (domesticCcy));
    SetForeignCcy(const_cast< ARM_Currency& > (foreignCcy));

    string cpnCcyName(cpnCcy.GetCcyName());

    string fundingCcyName(fundingCcy.GetCcyName());

    /// Set default keys for MDM data access
    ARM_StringVector keys(NbKeys);
    keys[YcKey]             = YC_KEY_NAME       + cpnCcyName;
    keys[OswModelKey]       = OSWMODEL_KEY_NAME + cpnCcyName;
    keys[CfModelKey]        = CFMODEL_KEY_NAME  + cpnCcyName;
    keys[MrsKey]            = MRS_KEY_NAME      + cpnCcyName;
 
	if( fundingCcyName == cpnCcyName )
    {
		keys.resize(NbKeys-1);
        // no basis, no forex but keep compatibility !
        keys[FundingKey]        = keys[YcKey];
        keys[BasisKey]          = keys[YcKey];
    }
    else
    {
        keys[FundingKey]        = YC_KEY_NAME       + fundingCcyName;
        keys[BasisKey]          = YC_BASIS_KEY_NAME + string(cpnCcy.GetCcyName())+"_"+ string(fundingCcy.GetCcyName());
        keys[ForexKey]          = FOREX_KEY_NAME + string(foreignCcy.GetCcyName()) + "_" + string(domesticCcy.GetCcyName());
    }

    SetKeys(keys);
 
    /// Set the coupon & funding model name alias i.e. keys to access to
    /// the right model for each GP keyword (1st argument)
    SetModelKeys();


    /// Create the Generic Security paid in coupon currency
    CreateAndSetDealDescriptionAndTimeIt(GetKeys()[itsCpnModelKey]);
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: Constructor
///	Returns: void
///	Action : builds the object (generic security only)
///			With Customized DateStrip
/////////////////////////////////////////////////////////////////
ARM_CRFCalculator::ARM_CRFCalculator(const ARM_DateStrip& fixLegSched,
					  const ARM_DateStrip& RFLegSched,
					  const ARM_DateStrip& ExerSched,
				const ARM_Date& asOfDate,
                const ARM_ReferenceValue& strike,
                int payRec,
                const ARM_Date& fixEndDate,
                int fixDayCount,
                int cpnDayCount,
                int cpnFreq,
                int cpnTiming,
                const string& cpnIndexTerm,
                int cpnIndexDayCount,
                const string& cpnResetCal,
                const string& cpnPayCal,
				int stubRule,
                int cpnResetGap,
                const ARM_ReferenceValue& leverage,
                const ARM_ReferenceValue& cpnMin,
                const ARM_ReferenceValue& cpnMax,
                const ARM_ReferenceValue& fundSpread,
                int fundFreq,
                int fundDayCount,
                const ARM_ReferenceValue& nominal,
                const ARM_ReferenceValue& exerFee,
                const ARM_Currency& cpnCcy,     
                const ARM_Currency& fundingCcy, 
                const ARM_Currency& basisCcy,   
                const ARM_Currency& domesticCcy,
                const ARM_Currency& foreignCcy, 
                const ARM_ReferenceValue& fundnominal,
				ARM_ModelType modelType)
:	ARM_GenCalculator(asOfDate),
    itsStrike(strike),
    itsPayRec(payRec),
    itsFixEndDate(fixEndDate),
    itsFixDayCount(fixDayCount),
    itsCpnDayCount(cpnDayCount),
    itsCpnFreq(cpnFreq),
    itsCpnTiming(cpnTiming),
    itsCpnIndexTerm(cpnIndexTerm),
    itsCpnIndexDayCount(cpnIndexDayCount),
    itsCpnResetCal(cpnResetCal),
    itsCpnPayCal(cpnResetCal),
	itsStubRule(stubRule),
    itsCpnResetGap(cpnResetGap),
    itsLeverage(leverage),
    itsCpnMin(cpnMin),
    itsCpnMax(cpnMax),
    itsFundSpread(fundSpread),
    itsFundFreq(fundFreq),
    itsFundDayCount(fundDayCount),
    itsNominal(nominal),
    itsFundNominal(fundnominal),
    itsExerGap(NOT_DEF_EXER_GAP),
    itsNbNonCall(NOT_DEF_EXER_GAP),
    itsExerFee(exerFee),
    itsOSWCalibFlag(ARM_SigmaCalibrationType::strikeEquivalent),
    itsMRSCalibType(ARM_MRSCalibrationType::Unknown),
    itsCapCalibFlag(true),
    itsFloorCalibFlag(false),
    itsProductToPrice(CRF_PRICE),
	itsCapFloorPF(NULL),
    itsFirstCallIdx(0),
	itsOWSSkewPF(NULL),
    itsModelType(modelType),
	itsMRSStrikeType(ARM_MRSStrikeCalibrationType::strikeEquivalent),
	itsSkewCalFlag(true),
	itsSkewReCalFlag(false),
	itsIsFrontier(0),
	itsExerCal(cpnResetCal),
	itsKeepOSWIdx(),
	itsvStrike()
{
    SetName(ARM_CALLREVFLOATER);

	itsStartDate = (ARM_Date) (ExerSched.GetFlowStartDates()->Elt(0));
	itsEndDate = (ARM_Date) (ExerSched.GetFlowEndDates()->Elt(ExerSched.GetFlowEndDates()->size()-1));

    /// Set the coupon/payment currency (inherited from ARM_Security)
    SetCurrencyUnit(const_cast< ARM_Currency* >(&cpnCcy));
    SetFundingCcy(const_cast< ARM_Currency& >(fundingCcy));
    SetBasisCcy(const_cast< ARM_Currency& >( basisCcy));
    SetDomesticCcy(const_cast< ARM_Currency& > (domesticCcy));
    SetForeignCcy(const_cast< ARM_Currency& > (foreignCcy));

    string cpnCcyName(cpnCcy.GetCcyName());

    string fundingCcyName(fundingCcy.GetCcyName());

    /// Set default keys for MDM data access
    ARM_StringVector keys(NbKeys);
    keys[YcKey]             = YC_KEY_NAME       + cpnCcyName;
    keys[OswModelKey]       = OSWMODEL_KEY_NAME + cpnCcyName;
    keys[CfModelKey]        = CFMODEL_KEY_NAME  + cpnCcyName;
    keys[MrsKey]            = MRS_KEY_NAME      + cpnCcyName;
 
	if( fundingCcyName == cpnCcyName )
    {
		keys.resize(NbKeys-1);
        // no basis, no forex but keep compatibility !
        keys[FundingKey]        = keys[YcKey];
        keys[BasisKey]          = keys[YcKey];
    }
    else
    {
        keys[FundingKey]        = YC_KEY_NAME       + fundingCcyName;
        keys[BasisKey]          = YC_BASIS_KEY_NAME + string(cpnCcy.GetCcyName())+"_"+ string(fundingCcy.GetCcyName());
        keys[ForexKey]          = FOREX_KEY_NAME + string(foreignCcy.GetCcyName()) + "_" + string(domesticCcy.GetCcyName());
    }

    SetKeys(keys);
 
    /// Set the coupon & funding model name alias i.e. keys to access to
    /// the right model for each GP keyword (1st argument)
    SetModelKeys();


    /// Create the Generic Security paid in coupon currency
	ARM_DateStripVector dateStrips;
	dateStrips.push_back((ARM_DateStrip*)(&fixLegSched));
	dateStrips.push_back((ARM_DateStrip*)(&RFLegSched));
	dateStrips.push_back((ARM_DateStrip*)(&ExerSched));
    CreateAndSetCustomDealDescriptionAndTimeIt(dateStrips,GetKeys()[itsCpnModelKey]);
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: Destructor
///	Returns: void
///	Action : destroys the object
/////////////////////////////////////////////////////////////////
ARM_CRFCalculator::~ARM_CRFCalculator()
{
    size_t i;
    for(i=0;i<itsStdSwaps.size();++i)
    {
        delete itsStdSwaps[i].first;
        itsStdSwaps[i].first=NULL;
    }

    for(i=0;i<itsFundStdSwaps.size();++i)
    {
        delete itsFundStdSwaps[i].first;
        itsFundStdSwaps[i].first=NULL;
    }

    for(i=0;i<itsIndexStdSwaps.size();++i)
    {
        delete itsIndexStdSwaps[i].first;
        itsIndexStdSwaps[i].first=NULL;
    }

	DeletePointorVector<ARM_VanillaArg>(itsVanillaArgVect); 

    delete itsCapFloorPF;
    itsCapFloorPF=NULL;
	delete 	itsOWSSkewPF;
	itsOWSSkewPF = NULL;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: Set???ToPrice() and Is???ToPrice()
///	Returns: void/boolean
///	Action : Set the flag to know which product to price.
///          Test the product currently able to be priced
/////////////////////////////////////////////////////////////////
void ARM_CRFCalculator::SetCRFToPrice()             {itsProductToPrice          =   CRF_PRICE;}
bool ARM_CRFCalculator::IsCRFToPrice() const        {return itsProductToPrice   ==  CRF_PRICE;}
void ARM_CRFCalculator::SetCapToPrice()             {itsProductToPrice          =   CAP_PRICE;}
bool ARM_CRFCalculator::IsCapToPrice() const        {return itsProductToPrice   ==  CAP_PRICE;}
void ARM_CRFCalculator::SetFloorToPrice()           {itsProductToPrice          =   FLOOR_PRICE;}
bool ARM_CRFCalculator::IsFloorToPrice() const      {return itsProductToPrice   ==  FLOOR_PRICE;}
void ARM_CRFCalculator::SetFundingToPrice()         {itsProductToPrice          =   FUNDING_PRICE;}
bool ARM_CRFCalculator::IsFundingToPrice() const    {return itsProductToPrice   ==  FUNDING_PRICE;}
void ARM_CRFCalculator::SetStdLegToPrice()          {itsProductToPrice          =   STDLEG_PRICE;}
bool ARM_CRFCalculator::IsStdLegToPrice() const     {return itsProductToPrice   ==  STDLEG_PRICE;}
void ARM_CRFCalculator::SetRFLegToPrice()           {itsProductToPrice          =   RFLEG_PRICE;}
bool ARM_CRFCalculator::IsRFLegToPrice() const      {return itsProductToPrice   ==  RFLEG_PRICE;}
void ARM_CRFCalculator::SetStdSwapToPrice()         {itsProductToPrice          =   STDSWAP_PRICE;}
bool ARM_CRFCalculator::IsStdSwapToPrice() const    {return itsProductToPrice   ==  STDSWAP_PRICE;}
void ARM_CRFCalculator::SetRFSwapToPrice()          {itsProductToPrice          =   RFSWAP_PRICE;}
bool ARM_CRFCalculator::IsRFSwapToPrice() const     {return itsProductToPrice   ==  RFSWAP_PRICE;}
void ARM_CRFCalculator::SetBermudaToPrice()         {itsProductToPrice          =   BERMUDA_PRICE;}
bool ARM_CRFCalculator::IsBermudaToPrice() const    {return itsProductToPrice   ==  BERMUDA_PRICE;}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: GetMRS & SetMRS
///	Returns: 
///	Action : get & set the MRS param for SUMMIT interface
/////////////////////////////////////////////////////////////////
const ARM_CurveModelParam& ARM_CRFCalculator::GetMRS() const
{
    return *(static_cast< ARM_CurveModelParam* >(GetMktDataManager()->GetData(GetKeys()[MrsKey])));
}

void ARM_CRFCalculator::SetMRS(ARM_ModelParam* mrsParam)
{
    if(!mrsParam || mrsParam->GetType() != ARM_ModelParamType::MeanReversion)
		    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : an MRS Param is expected !" );

    GetMktDataManager()->RegisterData(GetKeys()[MrsKey],static_cast< ARM_Object* >(mrsParam));
}

void ARM_CRFCalculator::SetMRS(ARM_ReferenceValue* mrsValues)
{
	if(mrsValues)
	{
		double	mrsValue = InterpolMeanRevParam(mrsValues);
		
		SetMRS(mrsValue);
	}
}

void ARM_CRFCalculator::SetMRS(double mrsValue)
{
	ARM_GP_Vector values(1, mrsValue);
	ARM_GP_Vector terms(1, 0.0);

	ARM_CurveModelParam * mrsParam = new ARM_CurveModelParam( ARM_ModelParamType::MeanReversion, mrsValue,"MRS");

    SetMRS(mrsParam);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: InitCRFForSummit
///	Returns: nothing
///	Action : Feed the MDM with input market objects
/////////////////////////////////////////////////////////////////
void ARM_CRFCalculator::InitCRFForSummit(ARM_ZeroCurve* zcCpn, 
										 ARM_VolCurve* swoptVC, 
										 ARM_VolCurve* capVC, 
										 ARM_VolLInterpol* capRho, 
										 ARM_VolLInterpol* capNu,
										 ARM_VolLInterpol* capBeta,
										 ARM_VolLInterpol* SwoptRho,
										 ARM_VolLInterpol* SwoptNu,
										 ARM_VolLInterpol* SwoptBeta,
										 ARM_ZeroCurve* zcFund, 
 										 ARM_ZeroCurve* zcCpnBasis, 
 										 ARM_ZeroCurve* zcFundBasis, 
										 double fxSpot,
                                         int UPDATE, // By default No (=0)
                                         long SkewReCalFlag,// By default(-1): IRRELEVANT!
                                         int SABRSigmaOrAlpha) // SABR Sigma or Alpha
                                                               // By default=1(Sigma)
{
    if ( SkewReCalFlag >= 0 )
    {
       if (SkewReCalFlag)
          itsSkewReCalFlag = true;
       else
          itsSkewReCalFlag = false;
    }
    // else let it as in the constructor
 
    /// Give a default value to MRS if not set.
    /// Set the STM calib flag
    ARM_CurveModelParam mr;

	vector <ARM_Object*> marketDatas;

    if (GetMktDataManager()->TestIfKeyMissing(GetKeys()[MrsKey]))    {
        ARM_GP_Vector values(1,MRS_DEFAULT_VALUE);
        ARM_GP_Vector times(1,0.0);
        mr = ARM_CurveModelParam(ARM_ModelParamType::MeanReversion,&values,&times,"MRS");
        itsMRSCalibType = ARM_MRSCalibrationType::stmfirstcolumn;
    }
    else
    {
       if ( UPDATE != 1 )
       {
          itsMRSCalibType = ARM_MRSCalibrationType::Unknown;
	      mr = GetMRS();
       }
       else
       {
          mr = GetMRS(); // itsShortTermCalibFlag already setted let it as is !!!
       }
    }

	marketDatas.push_back(zcCpn);

	ARM_BSModel* capBSmod = NULL;
	ARM_BSModel* SwoptBSmod = NULL;
	int SABR_Flag;
	
	if (capRho && capNu)
	{
		if (capBeta)
		{		
			SABR_Flag = capBeta->IsEqualToOne()? K_SABR_ARITH:K_SABR_IMPLNVOL;

			capBSmod = new ARM_BSSmiledModel(zcCpn->GetAsOfDate(), 
                                             0, 
                                             zcCpn, 
                                             zcCpn, 
                                             capVC, 
                                             K_YIELD, 
                                             capRho,
											 capNu, 
                                             SABR_Flag, 
                                             capBeta,
                                             0.5, // SABR Weight: Irrelevant
                                             SABRSigmaOrAlpha);
		}
		else
		{
			capBSmod = new ARM_BSSmiledModel(zcCpn->GetAsOfDate(), 
                                             0, 
                                             zcCpn, 
                                             zcCpn, 
                                             capVC, 
                                             K_YIELD, 
                                             capRho,
											 capNu, 
                                             K_SABR_ARITH,
                                             NULL,
                                             0.5, // SABR Weight: Irrelevant
                                             SABRSigmaOrAlpha);
		}
	}
	else
	{
		capBSmod = new ARM_BSModel(zcCpn->GetAsOfDate(), 0.0, zcCpn, zcCpn, capVC, K_YIELD);
	}

	if( SwoptRho && SwoptNu )
	{
		if (SwoptBeta)
		{
			SABR_Flag = SwoptBeta->IsEqualToOne()? K_SABR_ARITH : K_SABR_IMPLNVOL;
		
            SwoptBSmod = new ARM_BSSmiledModel(zcCpn->GetAsOfDate(), 
                                               0.0,
                                               zcCpn,
                                               zcCpn,
                                               swoptVC,
                                               K_YIELD,
										       SwoptRho,
                                               SwoptNu,
                                               SABR_Flag, 
                                               SwoptBeta,
                                               0.5, // SABR Weight: Irrelevant
                                               SABRSigmaOrAlpha);
		}
		else
		{
			SwoptBSmod = new ARM_BSSmiledModel(zcCpn->GetAsOfDate(), 
                                               0.0,
                                               zcCpn,
                                               zcCpn,
                                               swoptVC,
                                               K_YIELD,
										       SwoptRho,
                                               SwoptNu,
                                               K_SABR_ARITH,
                                               NULL,
                                               0.5, // SABR Weight: Irrelevant
                                               SABRSigmaOrAlpha);
		}
	}
	else
	{
		SwoptBSmod = new ARM_BSModel(zcCpn->GetAsOfDate(), 0.0, zcCpn, zcCpn, swoptVC, K_YIELD);
	}
										   

	marketDatas.push_back(SwoptBSmod);
	marketDatas.push_back(capBSmod);
	marketDatas.push_back(&mr);

	ARM_ZeroCurve* zcForDomBS = NULL;
	ARM_Forex* forex = NULL;

	// Update Mkt data when Basis CRF
	// Be careful, an order must be respected as in Summit Constructor
	if (IsBasis())
	{
		forex = new ARM_Forex(&GetForeignCcy(), &GetDomesticCcy(), fxSpot);

		///////*** APPEL FCT DEV PAR MAB !!! ***/////////
		// pour pouvoir tenir compte du basis de la jambe funding...
		// test pour savoir qui est USD

		ARM_CRV_TERMS psMatu;
		zcForDomBS = GenerateTwoCcyBSAdjusted(zcFund, zcFundBasis, 
                                              zcCpn, zcCpnBasis,
                                              0, psMatu);

		marketDatas.push_back(zcFund);
		marketDatas.push_back(zcForDomBS);
		marketDatas.push_back(forex);

    }

    if (UPDATE) // Just update correctly the calculator (in Hedge context)
       Update(marketDatas);
    else
       Init(marketDatas); // Creation of Market data manager

	if (zcForDomBS)
	   delete zcForDomBS;
	zcForDomBS = NULL;
	
	if (capBSmod)
	   delete capBSmod;
	capBSmod = NULL;

	if( SwoptBSmod )
		delete SwoptBSmod;
	SwoptBSmod = NULL;

	if (forex)
		delete forex;
	forex = NULL;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: SetModelKeys
///	Returns: void
///	Action : set the coupon, the funding & the basis reference
///          model names alias
/////////////////////////////////////////////////////////////////
void ARM_CRFCalculator::SetModelKeys()
{
    // Get the coupon & funding curve ccy
    string cpnCcyName( GetCurrencyUnit()->GetCcyName());
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
///	Class  : ARM_CRFCalculator
///	Routine: CheckData and CheckMktData
///	Returns: void
///	Action : check if CRF datas are consistent
/////////////////////////////////////////////////////////////////
void ARM_CRFCalculator::CheckData()
{
    if(itsNbNonCall < 0)
        itsNbNonCall = 0;

    /// For in-advance timing notification is necessary at cpn reset or before
    /// Here gaps are negative
    if(itsExerGap != NOT_DEF_EXER_GAP && itsCpnTiming == K_ADVANCE && itsExerGap > itsCpnResetGap)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : for in-advance timing, notification gap must be greater or equal than reset gap");

	CheckMktData();
}

void ARM_CRFCalculator::CheckMktData()
{
	/// MdM datas checking
	ARM_ZeroCurve* cpnCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));
    if(!cpnCurve)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Yield Curve for key =" + GetKeys()[YcKey] + " is expected in the Market Data Manager");
    string cpnCcy(cpnCurve->GetCurrencyUnit()->GetCcyName());


	ARM_ZeroCurve* fundingCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[FundingKey]));
    if(!fundingCurve)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Yield Curve for key =" + GetKeys()[FundingKey] + " is expected in the Market Data Manager");
    string fundingCcy(fundingCurve->GetCurrencyUnit()->GetCcyName());


	ARM_ZeroCurve* basisCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[BasisKey]));
    if(!basisCurve)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Yield Curve for key =" + GetKeys()[BasisKey] + " is expected in the Market Data Manager");
    string basisCcy(basisCurve->GetCurrencyUnit()->GetCcyName());

    if(cpnCcy != fundingCcy && basisCcy != cpnCcy && basisCcy != fundingCcy)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Basis Curve currency should be " + cpnCcy + " or " + fundingCcy);


    ARM_BSModel* oswBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswModelKey]) );
    if(!oswBSModel)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : swaption B&S model for key =" + GetKeys()[OswModelKey] + " is expected in the Market Data Manager");

    
    ARM_BSModel* cfBSModel = dynamic_cast< ARM_BSModel* >(GetMktDataManager()->GetData(GetKeys()[CfModelKey]) );
    if(!cfBSModel)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : cap/floor B&S model for key =" + GetKeys()[CfModelKey] + " is expected in the Market Data Manager");


	ARM_ModelParam* mrsParam = dynamic_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[MrsKey]));
    if(!mrsParam || mrsParam->GetType() != ARM_ModelParamType::MeanReversion)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : MRS Param for key =" + GetKeys()[MrsKey] + " is expected in the Market Data Manager");


    if(cpnCcy != fundingCcy)
    {
        if(GetMktDataManager()->TestIfKeyMissing(GetKeys()[ForexKey]))
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : a Forex is required for different coupon & funding currencies");

	    ARM_Forex* forex = dynamic_cast<ARM_Forex*>(GetMktDataManager()->GetData(GetKeys()[ForexKey]));
        if(!forex)
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Forex for key =" + GetKeys()[ForexKey] + " is expected in the Market Data Manager");

        string domesticCcy(forex->GetMainCurrency()->GetCcyName());
        if(domesticCcy != cpnCcy && domesticCcy != fundingCcy)
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Domestic currency should be " + cpnCcy + " or " + fundingCcy);

        string foreignCcy(forex->GetMoneyCurrency()->GetCcyName());
        if(foreignCcy != cpnCcy && foreignCcy != fundingCcy)
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Foreign currency should be " + cpnCcy + " or " + fundingCcy);
    }
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: ColumnNames
///	Returns: ARM_RowInfo
///	Action : create the column names of the deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_CRFCalculator::ColumnNames() const
{
    size_t colNamesSize = sizeof(CRFColNamesTable)/sizeof(CRFColNamesTable[0]);
    vector< string > colNamesVec(colNamesSize);
    vector< ARM_GP_VALUE_TYPE > colTypeVec(colNamesSize, ARM_STRING); 

    for(size_t i=0;i<colNamesSize; ++i)
        colNamesVec[i] = CRFColNamesTable[i];

    ARM_RowInfo rowInfo(colNamesVec,colTypeVec);

    return rowInfo;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: DatesStructure
///	Returns: ARM_DateStripCombiner
///	Action : create the list of all event dates of the CRF. The
///          DateStripCombiner merges event dates of each
///          legs
/////////////////////////////////////////////////////////////////
ARM_DateStripCombiner ARM_CRFCalculator::DatesStructure() const
{
    /// Get reset & payment calendars
    const char* resetCalendar   = itsCpnResetCal.c_str();
    const char* payCalendar     = itsCpnPayCal.c_str();
    const char* exerCalendar    = itsExerCal.c_str();

    int fwdRule		=   K_MOD_FOLLOWING;	// for forward dates
    int intRule		=   K_ADJUSTED;	// for interest dates
    
	int resetFreq   =   itsCpnFreq;
    int resetTiming =   itsCpnTiming;
	int payFreq     =   itsCpnFreq;
	int payGap      =   GETDEFAULTVALUE;
    int payTiming   =   K_ARREARS;

    ARM_Date actStartDate(itsStartDate);
    ARM_Date actFixEndDate(itsFixEndDate);


    /// Get standard reset gap
    int stdResetGap     = - GetCurrencyUnit()->GetSpotDays();
    int stdResetTiming  = K_ADVANCE;

    /// Build the fix coupon schedule with exercise dates
    ARM_DateStrip FixLegSched(actStartDate,actFixEndDate,resetFreq,itsFixDayCount,
            resetCalendar,
            fwdRule,intRule,itsStubRule,
            (itsExerGap==NOT_DEF_EXER_GAP ? stdResetGap : itsExerGap),
            payFreq,payGap,
            payCalendar,
            stdResetTiming,payTiming);

    /// Build the RF coupon schedule
    /// Always build with the cpn reset gap to get right forward start & end
    ARM_DateStrip RFLegSched(actFixEndDate,itsEndDate,resetFreq,itsCpnDayCount,
            resetCalendar,
            fwdRule,intRule,itsStubRule,
            itsCpnResetGap,
            payFreq,payGap,
            payCalendar,
            resetTiming,payTiming);

    /// Build the exercise schedule
    ARM_DateStrip ExerSched(actStartDate,itsEndDate,resetFreq,itsCpnDayCount,
            exerCalendar, // if exerCal has not been setted, it's the payCal
            fwdRule,intRule,itsStubRule,
            (itsExerGap==NOT_DEF_EXER_GAP ? stdResetGap : itsExerGap),
            payFreq,payGap,
            payCalendar,
            stdResetTiming,payTiming);

    ARM_DateStripVector SchedVect(NB_CRF_SCHED,NULL);
    SchedVect[FIX_CPN_SCHED]    = &FixLegSched;
    SchedVect[RF_CPN_SCHED]     = &RFLegSched;
    SchedVect[EXER_SCHED]       = &ExerSched;

    /// Synchronise "reset dates"
    ARM_GP_Vector* exerciseDates = ExerSched.GetResetDates();
    if(itsExerGap==NOT_DEF_EXER_GAP)
    {
        /// In the SUMMIT case, exercise dates are inputed
        /// then they replace all previously generated exercise dates

        /// Non callable dates will be insert and associated fees
        /// set to NON_CALL_FEE value
        //ReplaceExerDates(SchedVect[EXER_SCHED],SchedVect[RF_CPN_SCHED]);
		ReplaceExerDatesWithExerGap(SchedVect[EXER_SCHED]);

        /// Restore new exercise dates
        exerciseDates = ExerSched.GetResetDates();

        /// Replace reset dates of FixLegSched by the ExerSched ones to get
        /// only one event date per flow
        /// Loop forward because fix cpn & exercise legs are synchronised from start date
        ARM_GP_Vector* fixResetDates = FixLegSched.GetResetDates();
        if( fixResetDates && fixResetDates->size() && exerciseDates && exerciseDates->size() )
        {
            if(fixResetDates->size() > exerciseDates->size())
		        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : fix coupon number must be less than exercise number" );

            for(int i=0;i<fixResetDates->size();++i)
                (*fixResetDates)[i]=(*exerciseDates)[i];
        }
    }

    /// Replace reset dates of RFLegSched by the ExerSched ones to get
    /// only one event date per flow
    /// Loop backward because reverse cpn & exercise legs are synchronised from end date
    ARM_GP_Vector* cpnResetDates = RFLegSched.GetResetDates();
    if( cpnResetDates && cpnResetDates->size() && exerciseDates && exerciseDates->size() )
    {
        if(cpnResetDates->size() > exerciseDates->size())
		    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : reverse coupon number must be less than exercise number" );

        for(int i=cpnResetDates->size()-1,j=exerciseDates->size()-1;i>=0;--i,--j)
            (*cpnResetDates)[i]=(*exerciseDates)[j];
    }

    /// Merge schedules on "ResetDate"
    ARM_DateStripCombiner EventSchedule(SchedVect,"ResetDate");

    /// Initialise the memorisation of the first call line
    /// This line corresponds at less to the first event date
    /// posterior to the asOfDate
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	ARM_VectorPtr eventDates = EventSchedule.GetMergeData();
    itsFirstCallIdx = 0;
	while(itsFirstCallIdx < eventDates->size() && (*eventDates)[itsFirstCallIdx] <= asOfDate+K_NEW_DOUBLE_TOL)
        ++itsFirstCallIdx;

    return EventSchedule;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: DatesStructure
///	Returns: ARM_DateStripCombiner
///	Action : create the list of all event dates of the CRF. The
///          DateStripCombiner merges event dates of each
///          legs
/////////////////////////////////////////////////////////////////
ARM_DateStripCombiner ARM_CRFCalculator::CustomDatesStructure(const ARM_DateStripVector& dateStrips) const
{
	if (dateStrips.size() != NB_CRF_SCHED)
	{
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : invalid number of DateStrips" );
	}

    /// Get reset & payment calendars
    const char* resetCalendar   = itsCpnResetCal.c_str();
    const char* payCalendar     = itsCpnPayCal.c_str();
    const char* exerCalendar    = itsExerCal.c_str();

    int fwdRule		=   K_MOD_FOLLOWING;	// for forward dates
    int intRule		=   K_MOD_FOLLOWING;	// for interest dates
    
	int resetFreq   =   itsCpnFreq;
    int resetTiming =   itsCpnTiming;
	int payFreq     =   itsCpnFreq;
	int payGap      =   GETDEFAULTVALUE;
    int payTiming   =   K_ARREARS;

    ARM_Date actStartDate(itsStartDate);
    ARM_Date actFixEndDate(itsFixEndDate);


    /// Get standard reset gap
    int stdResetGap     = - GetCurrencyUnit()->GetSpotDays();
    int stdResetTiming  = K_ADVANCE;

    ARM_DateStrip FixLegSched( *dateStrips[FIX_CPN_SCHED]);

    ARM_DateStrip RFLegSched( *dateStrips[RF_CPN_SCHED]);

    ARM_DateStrip ExerSched( *dateStrips[EXER_SCHED]);

    ARM_DateStripVector SchedVect(NB_CRF_SCHED,NULL);
    SchedVect[FIX_CPN_SCHED]    = &FixLegSched;
    SchedVect[RF_CPN_SCHED]     = &RFLegSched;
    SchedVect[EXER_SCHED]       = &ExerSched;

    /// Synchronise "reset dates"
    ARM_GP_Vector* exerciseDates = ExerSched.GetResetDates();
    if(itsExerGap==NOT_DEF_EXER_GAP)
    {
        /// In the SUMMIT case, exercise dates are inputed
        /// then they replace all previously generated exercise dates

        /// Non callable dates will be insert and associated fees
        /// set to NON_CALL_FEE value
        //ReplaceExerDates(SchedVect[EXER_SCHED],SchedVect[RF_CPN_SCHED]);
		ReplaceExerDatesWithExerGap(SchedVect[EXER_SCHED]);

        /// Restore new exercise dates
        exerciseDates = ExerSched.GetResetDates();

        /// Replace reset dates of FixLegSched by the ExerSched ones to get
        /// only one event date per flow
        /// Loop forward because fix cpn & exercise legs are synchronised from start date
        ARM_GP_Vector* fixResetDates = FixLegSched.GetResetDates();
        if( fixResetDates && fixResetDates->size() && exerciseDates && exerciseDates->size() )
        {
            if(fixResetDates->size() > exerciseDates->size())
		        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : fix coupon number must be less than exercise number" );

            for(int i=0;i<fixResetDates->size();++i)
                (*fixResetDates)[i]=(*exerciseDates)[i];
        }
    }

    /// Replace reset dates of RFLegSched by the ExerSched ones to get
    /// only one event date per flow
    /// Loop backward because reverse cpn & exercise legs are synchronised from end date
    ARM_GP_Vector* cpnResetDates = RFLegSched.GetResetDates();
    if( cpnResetDates && cpnResetDates->size() && exerciseDates && exerciseDates->size() )
    {
        if(cpnResetDates->size() > exerciseDates->size())
		    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : reverse coupon number must be less than exercise number" );

        for(int i=cpnResetDates->size()-1,j=exerciseDates->size()-1;i>=0;--i,--j)
            (*cpnResetDates)[i]=(*exerciseDates)[j];
    }

    /// Merge schedules on "ResetDate"
    ARM_DateStripCombiner EventSchedule(SchedVect,"ResetDate");

    /// Initialise the memorisation of the first call line
    /// This line corresponds at less to the first event date
    /// posterior to the asOfDate
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	ARM_VectorPtr eventDates = EventSchedule.GetMergeData();
    itsFirstCallIdx = 0;
	while(itsFirstCallIdx < eventDates->size() && (*eventDates)[itsFirstCallIdx] <= asOfDate+K_NEW_DOUBLE_TOL)
        ++itsFirstCallIdx;

    return EventSchedule;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: ReplaceExerDates
///	Returns: void
///	Action : replace the reset dates of the date strip by
///          the exercise dates of the CRF and update fee table
///          to disable exercise at non relevant dates
/////////////////////////////////////////////////////////////////
void ARM_CRFCalculator::ReplaceExerDates(ARM_DateStrip* schedule,const ARM_DateStrip* cpnSchedule) const
{
    ARM_GP_Vector* resetDates = schedule->GetResetDates();
    if ( !resetDates || !resetDates->size() )
        return; // nothing to do !
	
	ARM_Vector* tmpExerDates  = itsExerFee.GetDiscreteDates();
    ARM_GP_Vector* exerDates = To_pARM_GP_Vector(tmpExerDates);
	CC_NS(std,auto_ptr)<ARM_GP_Vector> exerDatesTmp(exerDates);

    if(exerDates->size()==0)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : At least one exercise date is needed in the fees schedule" );

    ARM_GP_Vector newResetDates;
    ARM_GP_Vector nonCallResetDates;

    int i,resetIdx=0;

    /// Skip first exercise dates...
    int exerIdx=0;
    while(exerIdx<itsExerFee.GetSize() &&
        (*exerDates)[exerIdx] <= (*resetDates)[resetIdx] + K_NEW_DOUBLE_TOL)
        ++exerIdx;
    ///... to be just before or at the first reset date
    if( exerIdx > 0)
        --exerIdx;
    
    /// Use current exercise date and find associated reset date
    for(;exerIdx<itsExerFee.GetSize();++exerIdx)
    {
        /// Each intermediate reset date keeps its value but becomes non callable
        while(resetIdx < resetDates->size() &&
            (*resetDates)[resetIdx] < (*exerDates)[exerIdx] - K_NEW_DOUBLE_TOL)
        {
            nonCallResetDates.push_back((*resetDates)[resetIdx]);
            newResetDates.push_back((*resetDates)[resetIdx]);
            ++resetIdx;
        }

        if(resetIdx < resetDates->size())
        {
            /// Replace the current reset date by the exer date
            newResetDates.push_back((*exerDates)[exerIdx]);
            ++resetIdx;
        }
        else
            break;

    }

	delete exerDates;
	exerDates=NULL;

    /// For reset dates after the last exercise make them non callable
    /// and fit to the cpn reset date if anterior
    /// Loop backward to be synchronise with cpn schedule
    int firstNonCallIdx=resetIdx,cpnResetIdx;
    ARM_GP_Vector* cpnResetDates = cpnSchedule->GetResetDates();
    ARM_GP_Vector tmpResetDates;
    for(resetIdx=resetDates->size()-1,cpnResetIdx=cpnResetDates->size()-1;resetIdx>=firstNonCallIdx;--resetIdx,--cpnResetIdx)
    {
        if(cpnResetIdx>=0 && (*cpnResetDates)[cpnResetIdx] < (*resetDates)[resetIdx])
            tmpResetDates.push_back((*cpnResetDates)[cpnResetIdx]);
        else
            tmpResetDates.push_back((*resetDates)[resetIdx]);
    }
    /// Reverse the order to insert in increasing one !
    for(resetIdx=tmpResetDates.size()-1;resetIdx>=0;--resetIdx)
    {
        nonCallResetDates.push_back(tmpResetDates[resetIdx]);
        newResetDates.push_back(tmpResetDates[resetIdx]);
    }

    /// Replace reset dates in the target date strip
#ifdef __GP_STRICT_VALIDATION
    if(resetDates->size() != newResetDates.size())
		    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : consitency problem in merging exercise and reset schedules" );
#endif
    schedule->SetResetDates(&newResetDates);

    /// Update fees curve
    if(nonCallResetDates.size()>0)
    {
        double x,y;
        map<double,double> newExerFee;
        for(i=0;i<itsExerFee.GetSize();++i)
        {
            x=(*(itsExerFee.GetDiscreteDates()))[i];
            y=(*(itsExerFee.GetDiscreteValues()))[i];
            newExerFee.insert( pair<double,double>(x,y) );
        }
        y=NON_CALL_FEE;
        for(i=0;i<nonCallResetDates.size();++i)
        {
            x=nonCallResetDates[i];
            newExerFee.insert( pair<double,double>(x,y) );
        }
        map<double,double>::const_iterator it = newExerFee.begin();
        map<double,double>::const_iterator end = newExerFee.end();
        ARM_GP_Vector newDates(newExerFee.size());
        ARM_GP_Vector* newValues = new ARM_GP_Vector(newExerFee.size());
        for(i=0;it!=end;++it,++i)
        {
            newDates[i]=it->first;
            (*newValues)[i]=it->second;
        }
		ARM_Vector tmpNewDates	= To_ARM_Vector(newDates);
		ARM_Vector* tmpNewValues= To_pARM_Vector(newValues);
		delete newValues;
        itsExerFee.SetDiscreteDates( &tmpNewDates );  /// cloned
        itsExerFee.SetDiscreteValues(tmpNewValues); /// not cloned !!
    }

    /// Compute non call period from the fee table
    itsNbNonCall = 0;
    for(i=0;i<itsExerFee.size();++i)
    {
        if((*(itsExerFee.GetDiscreteValues()))[i] == NON_CALL_FEE)
            ++itsNbNonCall;
        else
            /// Stop  : we get the right number
            break;
    }
}

// similar to ReplaceExerDates but instead of getting reset and cpnSched,
// we compute exerGap from the first exercise date and the associated start date
void ARM_CRFCalculator::ReplaceExerDatesWithExerGap(ARM_DateStrip* schedule) const
{
    ARM_GP_Vector* resetDates = schedule->GetResetDates();
    if ( !resetDates || !resetDates->size() )
        return; // nothing to do !
	
	ARM_Vector* tmpExerDates  = itsExerFee.GetDiscreteDates();
    ARM_GP_Vector* exerDates = To_pARM_GP_Vector(tmpExerDates);
	CC_NS(std,auto_ptr)<ARM_GP_Vector> exerDatesPtr(exerDates);

    if(exerDates->size()==0)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : At least one exercise date is needed in the fees schedule" );

    ARM_GP_Vector newResetDates;
    ARM_GP_Vector nonCallResetDates;

    int i,resetIdx=0;

    /// Skip first exercise dates...
    int exerIdx=0;
    while(exerIdx<itsExerFee.GetSize() &&
        (*exerDates)[exerIdx] <= (*resetDates)[resetIdx] + K_NEW_DOUBLE_TOL)
        ++exerIdx;
    ///... to be just before or at the first reset date
    if( exerIdx > 0)
        --exerIdx;
    
	// calculate business days gap from exerDate and associated startDate
	int exerGap;
    ARM_GP_Vector* startDates = schedule->GetFlowStartDates();
	if (!startDates || (startDates->size() != resetDates->size()))
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : consitency problem in exercise sched between start and reset dates" );

	// get startDate associated with exerDate
	// not necessary the first, ex : 1st exerDate fits with 3rd reset period
	int tmpIdx = 0;
	while ( tmpIdx < resetDates->size() && 
			(*resetDates)[tmpIdx] < (*exerDates)[exerIdx] - K_NEW_DOUBLE_TOL )
		tmpIdx++;
		
	// calculate business days gap
	ARM_Date exerDate((*exerDates)[exerIdx]);
	ARM_Date startDate((*startDates)[tmpIdx]);
    char* exerCal = (char*)itsExerCal.c_str();
	if (tmpIdx < resetDates->size())
	{
		startDate.ChangeDate((*startDates)[tmpIdx]);
		exerGap = CountBusinessDays(startDate, exerDate, exerCal);
	}
	else
		exerGap = GETDEFAULTVALUE; // impossible, it would seems the 1st exerDate is after the last payDate ??

    /// Use current exercise date and find associated reset date
    for(;exerIdx<itsExerFee.GetSize();++exerIdx)
    {
        /// Each intermediate reset date keeps its value but becomes non callable
		// AV : we store (startDate - exerGap) instead of resetDate to be relevant with the first exercise date
        while(resetIdx < resetDates->size() &&
            (*resetDates)[resetIdx] < (*exerDates)[exerIdx] - K_NEW_DOUBLE_TOL)
        {
			exerDate.ChangeDate((*startDates)[resetIdx]);
			if (exerGap < 0)
				exerDate.PreviousBusinessDay(-exerGap, exerCal);
			else
				exerDate.NextBusinessDay(exerGap, exerCal);
            nonCallResetDates.push_back(exerDate.GetJulian());
            newResetDates.push_back(exerDate.GetJulian());
            ++resetIdx;
        }

        if(resetIdx < resetDates->size())
        {
            /// Replace the current reset date by the exer date
            newResetDates.push_back((*exerDates)[exerIdx]);
            ++resetIdx;
        }
        else
            break;
    }

    /// For reset dates after the last exercise make them non callable
	// AV : we store (startDate - exerGap) instead of resetDate to be relevant with the first exercise date
    for(;resetIdx < resetDates->size();resetIdx++)
    {
		exerDate.ChangeDate((*startDates)[resetIdx]);
		if (exerGap < 0)
			exerDate.PreviousBusinessDay(-exerGap, exerCal);
		else
			exerDate.NextBusinessDay(exerGap, exerCal);
        nonCallResetDates.push_back(exerDate.GetJulian());
        newResetDates.push_back(exerDate.GetJulian());
	}

    /// Replace reset dates in the target date strip
#ifdef __GP_STRICT_VALIDATION
    if(resetDates->size() != newResetDates.size())
		    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : consitency problem in merging exercise and reset schedules" );
#endif
    schedule->SetResetDates(&newResetDates);

    /// Update fees curve
    if(nonCallResetDates.size()>0)
    {
        double x,y;
        map<double,double> newExerFee;
        for(i=0;i<itsExerFee.GetSize();++i)
        {
            x=(*(itsExerFee.GetDiscreteDates()))[i];
            y=(*(itsExerFee.GetDiscreteValues()))[i];
            newExerFee.insert( pair<double,double>(x,y) );
        }
        y=NON_CALL_FEE;
        for(i=0;i<nonCallResetDates.size();++i)
        {
            x=nonCallResetDates[i];
            newExerFee.insert( pair<double,double>(x,y) );
        }
        map<double,double>::const_iterator it = newExerFee.begin();
        map<double,double>::const_iterator end = newExerFee.end();
        ARM_GP_Vector newDates(newExerFee.size());
        ARM_GP_Vector* newValues = new ARM_GP_Vector(newExerFee.size());
        for(i=0;it!=end;++it,++i)
        {
            newDates[i]=it->first;
            (*newValues)[i]=it->second;
        }
		ARM_Vector tmpNewDates	= To_ARM_Vector(newDates);
		ARM_Vector* tmpNewValues= To_pARM_Vector(newValues);
		delete newValues;
        itsExerFee.SetDiscreteDates( &tmpNewDates );  /// cloned
        itsExerFee.SetDiscreteValues(tmpNewValues); /// not cloned !!
    }

    /// Compute non call period from the fee table
    itsNbNonCall = 0;
    for(i=0;i<itsExerFee.size();++i)
    {
        if((*(itsExerFee.GetDiscreteValues()))[i] == NON_CALL_FEE)
            ++itsNbNonCall;
        else
            /// Stop  : we get the right number
            break;
    }
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: InitPriceableColumns
///	Returns: void
///	Action : initialise to 0 column to be able to sum
/////////////////////////////////////////////////////////////////
void ARM_CRFCalculator::InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const
{
    string zeroValue("0");

    /// Data auto-calibrated latter
    rowDescVec[AdjCapSpread] = zeroValue;
    rowTypeVec[AdjCapSpread] = ARM_DOUBLE;

    rowDescVec[AdjFloorSpread] = zeroValue;
    rowTypeVec[AdjFloorSpread] = ARM_DOUBLE;


    /// Set default 0 value for each column to be able to sum it
    rowDescVec[FundingVarFlow] = zeroValue;
    rowTypeVec[FundingVarFlow] = ARM_DOUBLE;

    rowDescVec[Funding] = zeroValue;
    rowTypeVec[Funding] = ARM_DOUBLE;

    rowDescVec[CpnIndexFlow] = zeroValue;
    rowTypeVec[CpnIndexFlow] = ARM_DOUBLE;

    rowDescVec[StdVarFlow] = zeroValue;
    rowTypeVec[StdVarFlow] = ARM_DOUBLE;

    rowDescVec[StdFixFlow] = zeroValue;
    rowTypeVec[StdFixFlow] = ARM_DOUBLE;

    rowDescVec[Caplet] = zeroValue;
    rowTypeVec[Caplet] = ARM_DOUBLE;

    rowDescVec[Floorlet] = zeroValue;
    rowTypeVec[Floorlet] = ARM_DOUBLE;

    rowDescVec[StdFlow] = zeroValue;
    rowTypeVec[StdFlow] = ARM_DOUBLE;

    rowDescVec[StdSwaplet] = zeroValue;
    rowTypeVec[StdSwaplet] = ARM_DOUBLE;

    rowDescVec[StdSwap] = zeroValue;
    rowTypeVec[StdSwap] = ARM_DOUBLE;

    rowDescVec[RFFlow] = zeroValue;
    rowTypeVec[RFFlow] = ARM_DOUBLE;

    rowDescVec[RFSwaplet] = zeroValue;
    rowTypeVec[RFSwaplet] = ARM_DOUBLE;

    rowDescVec[RFSwap] = zeroValue;
    rowTypeVec[RFSwap] = ARM_DOUBLE;

    rowDescVec[StdBermuda] = zeroValue;
    rowTypeVec[StdBermuda] = ARM_DOUBLE;

    rowDescVec[StdBermudaPrice] = zeroValue;
    rowTypeVec[StdBermudaPrice] = ARM_DOUBLE;

    rowDescVec[RFBermuda] = zeroValue;
    rowTypeVec[RFBermuda] = ARM_DOUBLE;

    rowDescVec[RFBermudaPrice] = zeroValue;
    rowTypeVec[RFBermudaPrice] = ARM_DOUBLE;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: GetLastValidDate
///	Returns: a line position in the deal description
///	Action : find the next available date in the deal description
///          for a input column name and a way to go further
/////////////////////////////////////////////////////////////////
int ARM_CRFCalculator::GetNextValidDateRow(const ARM_DealDescription& dealDescr,
                        CRFColAlias columnName,int prevRowIdx,int step,int stopIdx)
{
    int rowIdx=prevRowIdx;
    while( rowIdx != stopIdx && dealDescr.GetElemFormat(rowIdx,columnName) == ARM_MISSING_TYPE )
        rowIdx += step;

    return rowIdx;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: MiddleRows
///	Returns: ARM_RowInfo
///	Action : create a row of a deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_CRFCalculator::MiddleRows( size_t eventIdx, const ARM_DateStripCombiner& datesStructure ) const
{
    /// Here only exercise events are generated. At each date the coupon
    /// value is known with its analytical formula because
    /// keywords LIBOR & CAPLET were extended to handle the in arrears case
    /// For example of deal description with both coupon reset & bermuda notification
    /// dates, refer to previous version of 'crfcalculator.cpp' (before july 2004)
    size_t eventSize = datesStructure.GetDateStrip(0)->GetResetDates()->size();
    size_t descSize = sizeof(CRFColNamesTable)/sizeof(CRFColNamesTable[0]);

    vector< string > rowDescVec(descSize);
    vector< ARM_GP_VALUE_TYPE > rowTypeVec(descSize, ARM_MISSING_TYPE); 


    /// EventDate description
    double eventDate=(*(datesStructure.GetMergeData()))[eventIdx];
    CC_Ostringstream eventDateDesc;
    eventDateDesc << CC_NS(std,fixed) << eventDate;
    rowDescVec[EventDate] = eventDateDesc.str();
    rowTypeVec[EventDate] = ARM_DATE_TYPE;


    /// Convention conversion : K_QUATERLY => 3M for instance
    string fundMatFreq		= ARM_ArgConvReverse_MatFrequency.GetString(itsFundFreq);
    string fundDayCount		= ARM_ArgConvReverse_DayCount.GetString(itsFundDayCount);
    string cpnIndexDayCount	= ARM_ArgConvReverse_DayCount.GetString(itsCpnIndexDayCount);
    string cpnDayCount		= ARM_ArgConvReverse_DayCount.GetString(itsCpnDayCount);


    /// Set default 0 value for each column to be able to sum it
    InitPriceableColumns(rowDescVec,rowTypeVec);

    /// Get the model names for coupon & funding legs descritpion
    string cpnModelName = GetKeys()[itsCpnModelKey];
    string fundingModelName = GetKeys()[itsFundingModelKey];

    /// Get domestic & foreign model names in basis swap case for nominal conversion
    double spotAsOfValue = 1.0;
    string spotFutureValue("");
    if(fundingModelName != cpnModelName)
    {
        if( !GetMktDataManager()->TestIfKeyMissing(GetKeys()[ForexKey]) )
            spotAsOfValue = static_cast<ARM_Forex*>(GetMktDataManager()->GetData(GetKeys()[ForexKey]))->GetMarketPrice();
        if(string(GetFundingCcy().GetCcyName()) != string(GetDomesticCcy().GetCcyName()))
        {
            spotAsOfValue = 1.0/spotAsOfValue;
            spotFutureValue = "*SPOT(" + GetKeys()[ForexKey] + ")";
        }
        else
		{
            spotFutureValue = "/SPOT(" + GetKeys()[ForexKey] + ")";
		}
    }

    bool isLastExer = (eventIdx+1>=eventSize);

    /// Get the number of passed exercises
    size_t i,nbExerBefore=0;
    for(i=0;i<eventIdx;++i)
    {
        if((*(datesStructure.GetDateStrip(EXER_SCHED)->GetResetDates()))[i] != ARM_DateStripCombiner::DateStripCombiner_BlankData)
            ++nbExerBefore;
    }

    /// Define the right flow to describe the row
    string capletResetTiming("ADV");
    if(itsCpnTiming==K_ARREARS)
        capletResetTiming="ARR";

    /// Define the event type FIX or RF associated to the input row
    unsigned int eventType=0;
    if( (*(datesStructure.GetDateStrip(FIX_CPN_SCHED)->GetResetDates()))[eventIdx] != ARM_DateStripCombiner::DateStripCombiner_BlankData)
        eventType = FIX_CPN_EVENT;
    if( (*(datesStructure.GetDateStrip(RF_CPN_SCHED)->GetResetDates()))[eventIdx] != ARM_DateStripCombiner::DateStripCombiner_BlankData)
        eventType = RF_CPN_EVENT;

    string nextExerIdx("[i+1]");


    /// Fee description
    if(nbExerBefore < itsNbNonCall)
    {
        /// The call is not possible
        CC_Ostringstream feeDesc;
        feeDesc << CC_NS(std,fixed) << NON_CALL_FEE;
        rowDescVec[Fee] = feeDesc.str();
        itsFirstCallIdx = eventIdx + 1;
    }
    else
    {
        double exerFee;
        size_t exerFeeSize = itsExerFee.size();
        if(const_cast< ARM_ReferenceValue& >(itsExerFee).IsConstant() || exerFeeSize==0)
            exerFee = const_cast< ARM_ReferenceValue& >(itsExerFee).Interpolate(eventDate);
        else if(eventDate <= (*(const_cast< ARM_ReferenceValue& >(itsExerFee).GetDiscreteDates()))[exerFeeSize-1] + K_NEW_DOUBLE_TOL)
        {
            /// After the last notification date in the fees table,
            /// there isn't any possible exercise => fees are set to NON_CALL_FEE
            exerFee=const_cast< ARM_ReferenceValue& >(itsExerFee).Interpolate(eventDate);
        }
        else
            exerFee=NON_CALL_FEE;

        CC_Ostringstream feeDesc;
        feeDesc << exerFee;
        rowDescVec[Fee] = feeDesc.str();
    }
    rowTypeVec[Fee] = ARM_DOUBLE;


    /// Funding Flow description
    double  flowStartDate=(*(datesStructure.GetDateStrip(EXER_SCHED)->GetFlowStartDates()))[eventIdx];
    CC_Ostringstream fundStartDesc;
    fundStartDesc << CC_NS(std,fixed) << flowStartDate;
    rowDescVec[FundingStartDate] = fundStartDesc.str();
    rowTypeVec[FundingStartDate] = ARM_DATE_TYPE;

    double  flowEndDate=(*(datesStructure.GetDateStrip(EXER_SCHED)->GetFlowEndDates()))[eventIdx];
    CC_Ostringstream fundEndDesc;
    fundEndDesc << CC_NS(std,fixed) << flowEndDate;
    rowDescVec[FundingEndDate] = fundEndDesc.str();
    rowTypeVec[FundingEndDate] = ARM_DATE_TYPE;

    CC_Ostringstream fundSpreadDesc;
    fundSpreadDesc << CC_NS(std,fixed) << const_cast< ARM_ReferenceValue& >(itsFundSpread).Interpolate(flowStartDate);
    rowDescVec[FundingSpread] = fundSpreadDesc.str();
    rowTypeVec[FundingSpread] = ARM_DOUBLE;

    double flowPayDate=(*(datesStructure.GetDateStrip(EXER_SCHED)->GetPaymentDates()))[eventIdx];
    CC_Ostringstream fundNominalDesc;
    double nominalValue;
	// en cas de nominal cst size() renvoie 0, il faut tester directement sur itsDiscreteValues
    if(itsFundNominal.GetDiscreteValues() != NULL && itsFundNominal.GetDiscreteValues()->size() > 0)
    {
        nominalValue = const_cast< ARM_ReferenceValue& >(itsFundNominal).Interpolate(flowPayDate);
    }
    else
	{
        nominalValue = const_cast< ARM_ReferenceValue& >(itsNominal).Interpolate(flowPayDate);
        nominalValue *= spotAsOfValue;
	}

    fundNominalDesc << CC_NS(std,fixed) << CC_NS(std,setprecision)(2) <<nominalValue << spotFutureValue;
    rowDescVec[FundingNominal] = fundNominalDesc.str();
    rowTypeVec[FundingNominal] = ARM_STRING;


    CC_Ostringstream fundCommonDesc;
    fundCommonDesc << "(" << fundingModelName << "," << CRFColNamesTable[FundingStartDate] << "[i],";
    fundCommonDesc << CRFColNamesTable[FundingEndDate] << "[i],";
    fundCommonDesc << fundMatFreq << ",";
    fundCommonDesc << fundDayCount << ")";

    CC_Ostringstream fundAnnuityDesc;
    fundAnnuityDesc << "ANNUITY" << fundCommonDesc.str() << "*" << CRFColNamesTable[FundingNominal] << "[i]";
    rowDescVec[FundingAnnuity] = fundAnnuityDesc.str();
    rowTypeVec[FundingAnnuity] = ARM_STRING;

    CC_Ostringstream fundNominalExchangeDesc;
    fundNominalExchangeDesc << "(-DF(" << fundingModelName << "," << CRFColNamesTable[FundingStartDate] << "[i])";
    fundNominalExchangeDesc << "+DF(" << fundingModelName << "," << CRFColNamesTable[FundingEndDate] << "[i]))";
    fundNominalExchangeDesc << "*" << CRFColNamesTable[FundingNominal] << "[i]";
    fundNominalExchangeDesc << "+(DF(" << cpnModelName << "," << CRFColNamesTable[FundingStartDate] << "[i])";
    fundNominalExchangeDesc << "-DF(" << cpnModelName << "," << CRFColNamesTable[FundingEndDate] << "[i]))";
    fundNominalExchangeDesc << "*" << CRFColNamesTable[Nominal] << "[i]";
    rowDescVec[FundingNominalExchange] = fundNominalExchangeDesc.str();
    rowTypeVec[FundingNominalExchange] = ARM_STRING;

    CC_Ostringstream fundVarFlowDesc;
    fundVarFlowDesc << "SWAPRATE" << fundCommonDesc.str() << "*" << CRFColNamesTable[FundingAnnuity] << "[i]";
    fundVarFlowDesc << "+" << CRFColNamesTable[FundingNominalExchange] << "[i]";
    rowDescVec[FundingVarFlow] = fundVarFlowDesc.str();
    rowTypeVec[FundingVarFlow] = ARM_STRING;

    CC_Ostringstream fundSpreadFlowDesc;
    fundSpreadFlowDesc << CRFColNamesTable[FundingSpread] << "[i]*" << CRFColNamesTable[FundingAnnuity] << "[i]";
    rowDescVec[FundingSpreadFlow] = fundSpreadFlowDesc.str();
    rowTypeVec[FundingSpreadFlow] = ARM_STRING;


    CC_Ostringstream fundDesc;
    fundDesc << CRFColNamesTable[FundingVarFlow] << "[i]+" << CRFColNamesTable[FundingSpreadFlow]  << "[i]";
    rowDescVec[Funding] = fundDesc.str();
    rowTypeVec[Funding] = ARM_STRING;

        
    /// Underlying standard & reverse swaplets description (StdFlow or RFFlow vs Funding)
    CC_Ostringstream stdSwapletDesc;
    CC_Ostringstream rfSwapletDesc;
    if(eventIdx >= itsFirstCallIdx)
    {
        /// Std or RF Flow vs funding
        stdSwapletDesc << CRFColNamesTable[StdFlow] << "[i]";
        stdSwapletDesc << (itsPayRec==K_RCV ? "-" : "+") << CRFColNamesTable[Funding] << "[i]";
        rowDescVec[StdSwaplet] = stdSwapletDesc.str();
        rowTypeVec[StdSwaplet] = ARM_STRING;

        rfSwapletDesc << CRFColNamesTable[RFFlow] << "[i]";
        rfSwapletDesc << (itsPayRec==K_RCV ? "-" : "+") << CRFColNamesTable[Funding] << "[i]";
        rowDescVec[RFSwaplet] = rfSwapletDesc.str();
        rowTypeVec[RFSwaplet] = ARM_STRING;
    }


    /// Standard Swap & Option descriptions
    /// The swap is only described after the first call date
    CC_Ostringstream stdSwapDesc;
    if(eventIdx >= itsFirstCallIdx)
    {
        if(!isLastExer) 
            stdSwapDesc << "PV(" << CRFColNamesTable[StdSwap] << nextExerIdx << ")+";
        stdSwapDesc << CRFColNamesTable[StdSwaplet] << "[i]";
        rowDescVec[StdSwap] = stdSwapDesc.str();
        rowTypeVec[StdSwap] = ARM_STRING;
    }

    CC_Ostringstream stdBermudaDesc;
    stdBermudaDesc << "MAX(" << CRFColNamesTable[StdSwap] << "[i]-" << CRFColNamesTable[Fee] << "[i],";
    if(eventIdx == itsFirstCallIdx)
    {
        if(isLastExer)
            /// Last exercise => no PV()
            stdBermudaDesc << "0)";
        else
            /// Only one description to get the bermuda price
            stdBermudaDesc << "PV(" << CRFColNamesTable[StdBermuda] << nextExerIdx << "))";

        rowDescVec[StdBermudaPrice] = stdBermudaDesc.str();
        rowTypeVec[StdBermudaPrice] = ARM_STRING;
    }
    else if(eventIdx > itsFirstCallIdx)
    {
        if(isLastExer)
            /// Last exercise => no PV()
            stdBermudaDesc << "0)";
        else
            stdBermudaDesc << "PV(" << CRFColNamesTable[StdBermuda] << nextExerIdx << "))";
        rowDescVec[StdBermuda] = stdBermudaDesc.str();
        rowTypeVec[StdBermuda] = ARM_STRING;
    }


    /// RF Swap & Option description
    /// The swap is only described after the first call date
    CC_Ostringstream rfSwapDesc;
    if(eventIdx >= itsFirstCallIdx)
    {
        if(!isLastExer) 
            rfSwapDesc << "PV(" << CRFColNamesTable[RFSwap] << nextExerIdx << ")+";
        rfSwapDesc << CRFColNamesTable[RFSwaplet] << "[i]";
        rowDescVec[RFSwap] = rfSwapDesc.str();
        rowTypeVec[RFSwap] = ARM_STRING;
    }

	CC_Ostringstream finalDateStr;
    finalDateStr << CC_NS(std,fixed) << itsEndDate.GetJulian();
    rowDescVec[FinalDate] = finalDateStr.str();
    rowTypeVec[FinalDate] = ARM_DATE_TYPE;

	CC_Ostringstream SwapRateStr;
    SwapRateStr << "SwapRate(" << cpnModelName << "," << CRFColNamesTable[FundingStartDate] << "[i],";
    SwapRateStr << CRFColNamesTable[FinalDate] << "[i])";
	rowDescVec[SwapRate] = SwapRateStr.str();
	rowTypeVec[SwapRate] = ARM_STRING;


   
	CC_Ostringstream OptionToCallDesc;
    if(eventIdx == itsFirstCallIdx){
		OptionToCallDesc << CRFColNamesTable[RFBermuda]<< "[i]";
		rowDescVec[RFBermudaPrice] = OptionToCallDesc.str();
        rowTypeVec[RFBermudaPrice] = ARM_STRING;
	}

	CC_Ostringstream rfBermudaDesc;
	rfBermudaDesc << "MAX(" << CRFColNamesTable[RFSwap] << "[i]-" << CRFColNamesTable[Fee] << "[i],";
	if(isLastExer)
		rfBermudaDesc << "0)";
	else
		rfBermudaDesc << "PV(" << CRFColNamesTable[RFBermuda] << nextExerIdx << "))";
	rowDescVec[RFBermuda] = rfBermudaDesc.str();
	rowTypeVec[RFBermuda] = ARM_STRING;

	CC_Ostringstream rfFrontierDesc;
    rfFrontierDesc << "Frontier(" << CRFColNamesTable[RFSwap] << "[i]-" << CRFColNamesTable[Fee] << "[i],";
	if(isLastExer)
		rfFrontierDesc << "0,";
	else
		rfFrontierDesc << CRFColNamesTable[RFBermuda] << "[i+1],";

	rfFrontierDesc << CRFColNamesTable[SwapRate] << "[i])";
	rowDescVec[Frontier] = rfFrontierDesc.str();
	rowTypeVec[Frontier] = ARM_STRING;


    /// Flow Date descriptions
    int curSched = (eventType & FIX_CPN_EVENT) ? FIX_CPN_SCHED : RF_CPN_SCHED;
    flowStartDate = (*(datesStructure.GetDateStrip(curSched)->GetFlowStartDates()))[eventIdx];
    CC_Ostringstream startDateDesc;
    startDateDesc << CC_NS(std,fixed) << flowStartDate;
    rowDescVec[StartDate] = startDateDesc.str();
    rowTypeVec[StartDate] = ARM_DATE_TYPE;

    CC_Ostringstream endDateDesc;
    endDateDesc << CC_NS(std,fixed) << (*(datesStructure.GetDateStrip(curSched)->GetFlowEndDates()))[eventIdx];
    rowDescVec[EndDate] = endDateDesc.str();
    rowTypeVec[EndDate] = ARM_DATE_TYPE;

    CC_Ostringstream payDateDesc;
    double payDate=(*(datesStructure.GetDateStrip(curSched)->GetPaymentDates()))[eventIdx];
    payDateDesc << CC_NS(std,fixed) << payDate;
    rowDescVec[PayDate] = payDateDesc.str();
    rowTypeVec[PayDate] = ARM_DATE_TYPE;

    /// Payment DF description
    CC_Ostringstream dfPayDesc;
    dfPayDesc << "DF(" << cpnModelName << "," << CRFColNamesTable[PayDate] << "[i])";
    rowDescVec[DFPay] = dfPayDesc.str();
    rowTypeVec[DFPay] = ARM_STRING;

    /// Interest Period description 
    CC_Ostringstream itDesc;
    itDesc << CC_NS(std,fixed) << CC_NS(std,setprecision)(8) << (*(datesStructure.GetDateStrip(curSched)->GetInterestTerms()))[eventIdx];
    rowDescVec[IT] = itDesc.str();
    rowTypeVec[IT] = ARM_DOUBLE;

    /// Nominal description
    CC_Ostringstream nominalDesc;
    nominalDesc << CC_NS(std,fixed) << const_cast< ARM_ReferenceValue& >(itsNominal).Interpolate(payDate);
    rowDescVec[Nominal] = nominalDesc.str();
    rowTypeVec[Nominal] = ARM_DOUBLE;

    /// Strike description
    CC_Ostringstream strikeDesc;
    strikeDesc << CC_NS(std,fixed) << const_cast< ARM_ReferenceValue& >(itsStrike).Interpolate(flowStartDate);
    rowDescVec[Strike] = strikeDesc.str();
    rowTypeVec[Strike] = ARM_DOUBLE;
    double strike = atof(strikeDesc.str().c_str());

    /// Strike Flow description
    CC_Ostringstream strikeFlowDesc;
    strikeFlowDesc << CRFColNamesTable[Strike] << "[i]*";
    strikeFlowDesc << CRFColNamesTable[IT] << "[i]*" << CRFColNamesTable[Nominal] << "[i]*";
    strikeFlowDesc << CRFColNamesTable[DFPay] << "[i]";
    rowDescVec[StrikeFlow] = strikeFlowDesc.str();
    rowTypeVec[StrikeFlow] = ARM_STRING;

    /// Descriptions related to the fixed coupon event
    if(eventType & FIX_CPN_EVENT)
    {
        /// Fixed flow description (std flow is generated for bermuda swaption checking)
        CC_Ostringstream fixFlowDesc;
        fixFlowDesc << (itsPayRec==K_PAY ? "-" : " ") << CRFColNamesTable[StrikeFlow] << "[i]";
        rowDescVec[StdFlow] = fixFlowDesc.str();
        rowTypeVec[StdFlow] = ARM_STRING;

        CC_Ostringstream rfFlowDesc;
        rfFlowDesc << CRFColNamesTable[StdFlow] << "[i]";
        rowDescVec[RFFlow] = rfFlowDesc.str();
        rowTypeVec[RFFlow] = ARM_STRING;

    }// if fixed coupon event


    /// Descriptions related to the RF coupon event
    else if(eventType & RF_CPN_EVENT)
    {
        /// Cpn Index Date descriptions
        CC_Ostringstream indexStartDateDesc;
        indexStartDateDesc << CC_NS(std,fixed) << (*(datesStructure.GetDateStrip(RF_CPN_SCHED)->GetFwdRateStartDates()))[eventIdx];
        rowDescVec[IndexStartDate] = indexStartDateDesc.str();
        rowTypeVec[IndexStartDate] = ARM_DATE_TYPE;

        /// Get the standard gap of the cpn index
        int stdCpnIndexResetGap = -GetCurrencyUnit()->GetSpotDays();

        /// Cpn Index description
        CC_Ostringstream cpnIndexDesc;
        cpnIndexDesc << "LIBOR(" << cpnModelName << "," << CRFColNamesTable[IndexStartDate] << "[i],";
        cpnIndexDesc << itsCpnIndexTerm << "," << cpnIndexDayCount << ",";
        cpnIndexDesc << stdCpnIndexResetGap << "," << CRFColNamesTable[PayDate] << "[i])";
        rowDescVec[CpnIndex] = cpnIndexDesc.str();
        rowTypeVec[CpnIndex] = ARM_STRING;

        /// Coupon Min, Max & Leverage descriptions
        CC_Ostringstream cpnMinDesc;
        cpnMinDesc << CC_NS(std,fixed) << const_cast< ARM_ReferenceValue& >(itsCpnMin).Interpolate(flowStartDate);
        rowDescVec[CpnMin] = cpnMinDesc.str();
        rowTypeVec[CpnMin] = ARM_DOUBLE;
        double cpnMin = atof(cpnMinDesc.str().c_str());

        CC_Ostringstream cpnMaxDesc;
        cpnMaxDesc << CC_NS(std,fixed) << const_cast< ARM_ReferenceValue& >(itsCpnMax).Interpolate(flowStartDate);
        rowDescVec[CpnMax] = cpnMaxDesc.str();
        rowTypeVec[CpnMax] = ARM_DOUBLE;
        double cpnMax = atof(cpnMaxDesc.str().c_str());

        if(cpnMin > cpnMax)
        {
		    CC_Ostringstream os;
		    os << ARM_USERNAME << " : CpnMin > CpnMax EventDate = " << ARM_Date(eventDate).toString();
		    throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
        }

        CC_Ostringstream leverageDesc;
        leverageDesc << CC_NS(std,fixed) << const_cast< ARM_ReferenceValue& >(itsLeverage).Interpolate(flowStartDate);
        rowDescVec[Leverage] = leverageDesc.str();
        rowTypeVec[Leverage] = ARM_DOUBLE;
        double leverage=atof(leverageDesc.str().c_str());

        if(fabs(leverage) > K_NEW_DOUBLE_TOL)
        {
			if(leverage < -K_NEW_DOUBLE_TOL && itsOSWCalibFlag == ARM_SigmaCalibrationType::strikeEquivalent)
			{
				itsOSWCalibFlag = ARM_SigmaCalibrationType::StrikeDeltaApproxi;
			/*	CC_Ostringstream os;
				os << ARM_USERNAME << " : When leverage is negative the " << ARM_ArgConvReverse_SigmaCalibType.GetString( itsOSWCalibFlag );
				os << " type is not valid to calibrate sigma, please try DELTAPPROXI or ATM";
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );*/
			}

			string CAP   = leverage> K_NEW_DOUBLE_TOL ? "CAP":"FLOOR";
			string FLOOR = leverage> K_NEW_DOUBLE_TOL ? "FLOOR":"CAP";
			string sign   = leverage> K_NEW_DOUBLE_TOL ? " ":"-";
			/// Coupon index flow description
			CC_Ostringstream cpnIndexFlowDesc;
			cpnIndexFlowDesc << CRFColNamesTable[Leverage] << "[i]*" << CRFColNamesTable[CpnIndex] << "[i]*";
			cpnIndexFlowDesc << CRFColNamesTable[IT] << "[i]*" << CRFColNamesTable[DFPay] << "[i]*";
			cpnIndexFlowDesc << CRFColNamesTable[Nominal] << "[i]";
			rowDescVec[CpnIndexFlow] = cpnIndexFlowDesc.str();
			rowTypeVec[CpnIndexFlow] = ARM_STRING;
        
            if(cpnMin > strike)
            {
		        CC_Ostringstream os;
		        os << ARM_USERNAME << " : CpnMin > Strike at EventDate = " << ARM_Date(eventDate).toString();
		        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
            }

            /// Caplet description with analytical formula
            CC_Ostringstream capletDesc;
			double strikeMin = (strike - cpnMin)/leverage;
			if(strikeMin > K_NEW_DOUBLE_TOL){
				capletDesc << sign <<CRFColNamesTable[Leverage] << "[i]*CAPLET(" << cpnModelName << ",";
				capletDesc << CRFColNamesTable[StartDate] << "[i],";
				capletDesc << CRFColNamesTable[EndDate] << "[i],";			
				leverage > K_NEW_DOUBLE_TOL ? capletDesc << CRFColNamesTable[AdjCapSpread] << "[i]+(": capletDesc << CRFColNamesTable[AdjFloorSpread] << "[i]+(";
				capletDesc << CRFColNamesTable[Strike] << "[i]-";
				capletDesc << CRFColNamesTable[CpnMin] << "[i])/" << CRFColNamesTable[Leverage] << "[i],";
				capletDesc << CAP << "," << cpnDayCount << "," << itsCpnResetGap << ",PayDate[i],";
				capletDesc << CRFColNamesTable[Nominal] << "[i]," << capletResetTiming << ",";
				capletDesc << itsCpnIndexTerm << "," << cpnIndexDayCount << ")";
			}
			else if(CAP =="CAP")
			{
				capletDesc << "SWAP(" << cpnModelName << "," << CRFColNamesTable[StartDate] << "[i], "; 
				capletDesc << CRFColNamesTable[EndDate] << "[i], ";
				leverage > K_NEW_DOUBLE_TOL ? capletDesc << CRFColNamesTable[AdjCapSpread] << "[i]+(": capletDesc << CRFColNamesTable[AdjFloorSpread] << "[i]+(";
				capletDesc << CRFColNamesTable[Strike] << "[i]-";
				capletDesc << CRFColNamesTable[CpnMin] << "[i])/" << CRFColNamesTable[Leverage] << "[i],REC,,,";
				capletDesc << itsCpnIndexTerm << "," << cpnIndexDayCount <<",," ;
				capletDesc << CRFColNamesTable[Nominal] << "[i])";
			}
			else
			{
				capletDesc << 0.0;
			}
			rowDescVec[Caplet] = capletDesc.str();
			rowTypeVec[Caplet] = ARM_STRING;


            /// Floorlet description with analytical formula
            CC_Ostringstream floorletDesc;
			double strikeMax = (strike - cpnMax)/leverage;
			if(strikeMax > K_NEW_DOUBLE_TOL){
				floorletDesc << sign << CRFColNamesTable[Leverage] << "[i]*CAPLET(" << cpnModelName << ",";
				floorletDesc << CRFColNamesTable[StartDate] << "[i],";
				floorletDesc << CRFColNamesTable[EndDate] << "[i],";
				leverage > K_NEW_DOUBLE_TOL ? floorletDesc << CRFColNamesTable[AdjFloorSpread] << "[i]+(": floorletDesc << CRFColNamesTable[AdjCapSpread] << "[i]+(";
				floorletDesc << CRFColNamesTable[Strike] << "[i]-";
				floorletDesc << CRFColNamesTable[CpnMax] << "[i])/" << CRFColNamesTable[Leverage] << "[i],";
				floorletDesc << FLOOR << "," << cpnDayCount << "," << itsCpnResetGap << ",PayDate[i],";
				floorletDesc << CRFColNamesTable[Nominal] << "[i]," << capletResetTiming << ",";
				floorletDesc << itsCpnIndexTerm << "," << cpnIndexDayCount << ")";
				
			}
			else if(FLOOR =="FLOOR")
				{
					floorletDesc << 0.0;
				}
			else
			{
				floorletDesc << "SWAP(" << cpnModelName << "," << CRFColNamesTable[StartDate] << "[i], "; 
				floorletDesc << CRFColNamesTable[EndDate] << "[i], ";
				leverage > K_NEW_DOUBLE_TOL ? floorletDesc << CRFColNamesTable[AdjFloorSpread] << "[i]+(": floorletDesc << CRFColNamesTable[AdjCapSpread] << "[i]+(";
				floorletDesc << CRFColNamesTable[Strike] << "[i]-";
				floorletDesc << CRFColNamesTable[CpnMin] << "[i])/" << CRFColNamesTable[Leverage] << "[i], REC,,,";
				floorletDesc << itsCpnIndexTerm << "," << cpnIndexDayCount <<",," ;
				floorletDesc << CRFColNamesTable[Nominal] << "[i])";
			}
			rowDescVec[Floorlet] = floorletDesc.str();
			rowTypeVec[Floorlet] = ARM_STRING;
        }
        else
        {
		    CC_Ostringstream os;
		    os << ARM_USERNAME << " : Leverage null at EventDate = " << ARM_Date(eventDate).toString();
		    throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
        }

        /// Standard Flow description
        CC_Ostringstream stdFlowDesc;
        if(itsPayRec==K_RCV)
            stdFlowDesc << CRFColNamesTable[StrikeFlow] << "[i]-" << CRFColNamesTable[CpnIndexFlow] << "[i]";
        else
            stdFlowDesc << CRFColNamesTable[CpnIndexFlow] << "[i]-" << CRFColNamesTable[StrikeFlow] << "[i]";

        rowDescVec[StdFlow] = stdFlowDesc.str();
        rowTypeVec[StdFlow] = ARM_STRING;

        /// RF Flow description
        CC_Ostringstream rfFlowDesc;
        rfFlowDesc << CRFColNamesTable[StdFlow] << "[i]";
        if(itsPayRec==K_RCV)
            rfFlowDesc << "+" << CRFColNamesTable[Caplet] << "[i]-" << CRFColNamesTable[Floorlet] << "[i]";
        else
            rfFlowDesc << "-" << CRFColNamesTable[Caplet] << "[i]+" << CRFColNamesTable[Floorlet] << "[i]";
        rowDescVec[RFFlow] = rfFlowDesc.str();
        rowTypeVec[RFFlow] = ARM_STRING;

    }// if RF coupon event


    /// Standard var (cpn index & funding) & fixed (strike & funding spread)
    /// flows for calibration swaptions
    CC_Ostringstream stdFixFlowDesc;
    CC_Ostringstream stdVarFlowDesc;
    if(itsPayRec==K_PAY)
    {
        stdVarFlowDesc << CRFColNamesTable[FundingVarFlow] << "[i]+" << CRFColNamesTable[CpnIndexFlow] << "[i]";
        stdFixFlowDesc << CRFColNamesTable[FundingSpreadFlow] << "[i]-" << CRFColNamesTable[StrikeFlow] << "[i]";
    }
    else
    {
        stdVarFlowDesc << "-(" << CRFColNamesTable[FundingVarFlow] << "[i]+" << CRFColNamesTable[CpnIndexFlow] << "[i])";
        stdFixFlowDesc << CRFColNamesTable[StrikeFlow] << "[i]-" << CRFColNamesTable[FundingSpreadFlow] << "[i]";
    }
    rowDescVec[StdVarFlow] = stdVarFlowDesc.str();
    rowTypeVec[StdVarFlow] = ARM_STRING;
    rowDescVec[StdFixFlow] = stdFixFlowDesc.str();
    rowTypeVec[StdFixFlow] = ARM_STRING;

    return ARM_RowInfo(rowDescVec,rowTypeVec);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: CreateAndSetModel
///	Returns: void
///	Action : creates the model
/////////////////////////////////////////////////////////////////
void ARM_CRFCalculator::CreateAndSetModel()
{
	/// Get asOfDate
	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();
    /// Build the default stochastic model of the calculator : H&W 1F
	ARM_ZeroCurve* pCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));

	/// Creates default values for volatility & mean reversion
    /// (calibration will be called at pricing time) and set them
    /// in the reference model
	ARM_CurveModelParam volParam( ARM_ModelParamType::Volatility,SIGMA_DEFAULT_VALUE ,"SIGMA" );
	/// Get constant Mean Reversion from Curve Model Param.
	ARM_CurveModelParam* mrsParam = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(GetKeys()[MrsKey]));
	//double  flowStartDate=(*(datesStructure.GetDateStrip(EXER_SCHED)->GetFlowStartDates()))[itsFirstCallIdx];
	double tenorLag = itsStartDate > asOfDate ? itsEndDate - itsStartDate : itsEndDate - asOfDate;
	double mrsValue = mrsParam->GetValue(tenorLag);

	ARM_CurveModelParam CstmrsParam( ARM_ModelParamType::MeanReversion, mrsValue,"MRS");

	SetMRS(&CstmrsParam);
	ARM_ModelParamVector paramVector(2);
	paramVector[0] = &volParam;
	paramVector[1] = &CstmrsParam;

    /// Witch Model have we use to price
    ARM_PricingModelPtr refModel;
	switch(itsModelType)
	{
	case ARM_PricingModelType::HWM1F:
		refModel = ARM_PricingModelPtr(ARM_HWModelFactory.Instance()->CreateHWModel(CreateClonedPtr( pCurve ),paramVector));
		break;
	case ARM_PricingModelType::QGM1F:
        {
            ARM_CurveModelParam skewParam( ARM_ModelParamType::Skew,SKEW_DEFAULT_VALUE ,"SKEW" );
            paramVector.push_back(&skewParam);
            refModel = ARM_PricingModelPtr(static_cast< ARM_PricingModel* >(new ARM_QGM1F(CreateClonedPtr( pCurve ),paramVector)));
        }
        break;
	default:
		ARM_THROW( ERR_INVALID_ARGUMENT, "Only HWM / QGM are avaliable to price CRF" );
	
	}
       /// Build the model with or without basis effect
    /// build MultiAssets
    ARM_StringVector names (1);
	vector<ARM_PricingModelPtr> models (1);
	ARM_StringVectorVector depends(1);

    /// Key names
    names[myRefModel] = GetKeys()[YcKey];

    /// models
    models[myRefModel] = refModel;

    /// basis
    if(IsBasis())
    {
        /// Push back model pricing names
        names.resize(NbModels);
        names[myIrMarginModel] = GetKeys()[FundingKey];
        names[myBasisMarginModel] = GetKeys()[BasisKey];
        names[myForexModel] = GetKeys()[ForexKey];

        /// built modep pricing for multi-currency pricing models
        
        models.resize(NbModels);
        depends.resize(NbModels);
        /// Build the IR forward margin model
	    ARM_ZeroCurve* fundingCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[FundingKey]));
        models[myIrMarginModel] = ARM_PricingModelPtr(static_cast< ARM_PricingModel* >(new ARM_ForwardMarginIR( CreateClonedPtr(fundingCurve))));
        depends[myIrMarginModel] = ARM_StringVector (1,names[myRefModel]);

        /// Build the Basis forward margin model
	    ARM_ZeroCurve* basisCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[BasisKey]));
        if(itsBasisRefModelKey == YcKey)
        {
            models[myBasisMarginModel] = ARM_PricingModelPtr(static_cast< ARM_PricingModel* >(new ARM_ForwardMarginBasis( CreateClonedPtr(basisCurve)))); 
            depends[myBasisMarginModel] = ARM_StringVector (1,names[myRefModel]);
        }
        else{
            models[myBasisMarginModel] = ARM_PricingModelPtr(static_cast< ARM_PricingModel* >(new ARM_ForwardMarginBasis( CreateClonedPtr(basisCurve)))); 
            depends[myBasisMarginModel] = ARM_StringVector (1,names[myIrMarginModel]);
        }


        /// Build the Forward Forex model
	    ARM_Forex* forex = dynamic_cast<ARM_Forex*>(GetMktDataManager()->GetData(GetKeys()[ForexKey]));
        string domesticCcy(GetDomesticCcy().GetCcyName());
        string foreignCcy(GetForeignCcy().GetCcyName());
        string fundingCcy(GetFundingCcy().GetCcyName());
        string basisCcy(GetBasisCcy().GetCcyName());
        depends[myForexModel].resize(2);
        if(domesticCcy == basisCcy) 
        {
            models[myForexModel] = (foreignCcy == fundingCcy) ? ARM_PricingModelPtr( static_cast< ARM_PricingModel* >(
                                                     new ARM_ForwardForex(*forex,CreateClonedPtr( fundingCurve ),CreateClonedPtr( basisCurve)))) 
                                                   : ARM_PricingModelPtr( static_cast< ARM_PricingModel* >(
                                                     new ARM_ForwardForex(*forex,refModel->GetZeroCurve(),CreateClonedPtr( basisCurve))));
            depends[myForexModel][0] = names[myBasisMarginModel];
            depends[myForexModel][1] = (foreignCcy == fundingCcy) ? names[myIrMarginModel]: names[myRefModel];
        }
        else
        {
            models[myForexModel] = (domesticCcy == fundingCcy) ? ARM_PricingModelPtr( static_cast< ARM_PricingModel* >(
                                                     new ARM_ForwardForex(*forex,CreateClonedPtr( fundingCurve ),CreateClonedPtr( basisCurve))))
                                                   : ARM_PricingModelPtr( static_cast< ARM_PricingModel* >(
                                                     new ARM_ForwardForex(*forex,refModel->GetZeroCurve(),CreateClonedPtr( basisCurve))));
            depends[myForexModel][0] = (domesticCcy == fundingCcy) ? names[myIrMarginModel]: names[myRefModel];
            depends[myForexModel][1] = names[myBasisMarginModel];
        }
    }

	ARM_ModelNameMap modelMap (names, models, depends);
		
	///	modelMap & correls are cloned in multi assets model
	ARM_PricingModelPtr model (new ARM_MultiAssetsModel ( &modelMap) );

    /// Create a TreeMethod with a default step number per year and set it
    const ARM_DealDescription dealDescr = GetGenSecurity()->GetDealDescription();
    double lastEventTime = atof(dealDescr.GetElem(dealDescr.GetRowsNb()-1,EventDate).c_str()) - asOfDate.GetJulian();
    int nbSteps=static_cast<int>(floor(TREE_NBSTEPS_PER_YEAR*lastEventTime/K_YEAR_LEN));

	int schedulerType=ARM_SchedulerBase::ConstantVarianceMeanReverting;
	ARM_GP_Vector schedulerDatas(3);
	schedulerDatas[0] = nbSteps;
	schedulerDatas[1] = 1;
	schedulerDatas[2] = 1.0e-3;
	int samplerType=ARM_SamplerBase::MeanReverting;
	ARM_GP_Vector samplerDatas(1,1.0e-3);
	int truncatorType=ARM_TruncatorBase::StandardDeviation;
	ARM_GP_Vector truncatorDatas(1,5.0);
	int reconnectorType=ARM_ReconnectorBase::Mean;
	int smootherType=ARM_SmootherBase::DoNothing;
	ARM_TreeBase* tree = ARM_TreeFactory.Instance()->CreateTreeND(1,schedulerType,schedulerDatas,
		samplerType,samplerDatas,truncatorType,truncatorDatas,false,reconnectorType,smootherType);
	model->SetNumMethod( ARM_NumMethodPtr( tree ) );

    /// Create a Numeraire and set it
    ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::TerminalZc ) );
    model->SetNumeraire(numeraire);


	/// Set the model
	SetPricingModel(model);
}
/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: PreCalibrateSkew
///	Returns: ARM_CalibMethod
///	Action : to precalibrate  QG Skew
/////////////////////////////////////////////////////////////////
void ARM_CRFCalculator::PreCalibrateSkew()
{
	ARM_ZeroCurve* pCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));
	ARM_CurveModelParam* mrsParam = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(GetKeys()[MrsKey]));
	ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();
	ARM_BSModel* oswBSModel = GetOSWBSModel();
	size_t sizepf = GetOSWPortfolio()->size();
	ARM_BSModel flatbsModel(pCurve,NULL);
	double defaultATMK = 2.50;
	double initVol = 0.0;
	vector< ARM_Security* > swaptionList(1);
	vector< double> prices(1);
	vector< double> vegas(1);
	vector< double> weight(1,1.0);
	ARM_CurveModelParam  modelParamMrs(ARM_ModelParamType::MeanReversion, 0.0,"MRS");
	ARM_CurveModelParam  modelParamVol(ARM_ModelParamType::Volatility, initVol,"SIGMA",SIGMA_MIN_INIT_VALUE,SIGMA_MAX_INIT_VALUE);
	ARM_ModelParamVector volmodelparamVect(1,&modelParamVol);
	ARM_CalibMethod volCalibMethod(ARM_StdPortfolioPtr(NULL),volmodelparamVect,ARM_CalibMethodType::Bootstrap1D,MAXITER);

	ARM_CurveModelParam  modelParamSkew(ARM_ModelParamType::Skew, SKEW_DEFAULT_VALUE,"SKEW",SKEW_LOWER_BOUND,SKEW_UPPER_BOUND);
	ARM_ModelParamVector skewmodelparamVect(1,&modelParamSkew);
	ARM_ModelFitterDes modelfitter(ARM_ModelFitterSolverType::NagSolver,MAXITER,SKEW_CALIB_PRECISION,SKEW_CALIB_PRECISION);
	ARM_CalibMethod skewCalibMethod(ARM_StdPortfolioPtr(NULL),skewmodelparamVect,ARM_CalibMethodType::Bootstrap1D,&modelfitter,
		ARM_CalibrationTarget::PriceTarget,&volCalibMethod,NULL,true);

	ARM_ModelParamVector paramVector(3);
	paramVector[0] = &modelParamVol;
	paramVector[1] = &modelParamSkew;
	paramVector[2] = &modelParamMrs;
	ARM_QGM1F QGmodel(CreateClonedPtr( pCurve ),paramVector);
	
	/// creates the default numeraire for the time being
	ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::TerminalZc ) );
	QGmodel.SetNumeraire( numeraire );

	ARM_GP_Vector SkewTimes(sizepf),
		 VolTimes(sizepf),	
		 SkewValues(sizepf),
		 VolValues(sizepf),
		 weights(sizepf);
	double Initweight = 300.0;
	double amortweight = 0.8;
	vector< ARM_Security* > swaptionpfList(sizepf);
	vector< double> pfprices(sizepf);
	vector< double> pfvegas(sizepf);
	vector< double> pfweight(sizepf,1.0);
	double coefUp = 1.5;
	double coefDown = 0.5;
	
	for(int i=0; i<sizepf; i++)
	{
        ARM_Swaption* swaption= (ARM_Swaption*) GetOSWPortfolio()->GetAsset(i)->Clone();
        swaption->SetModel(oswBSModel);
        double SwapRate = swaption->ComputeBSSpot(swaption->GetExpiryDate());
        swaption->SetModel(NULL);
        double Keq = swaption->GetStrike();
		double skewStrike,sigmaStrike;

		if(itsIsFrontier){
			skewStrike = CC_Max(itsvStrike[i],Keq);
			sigmaStrike = CC_Min(itsvStrike[i],Keq);		
		}
		else{
			skewStrike = Keq > coefUp*SwapRate ? 0.75*SwapRate : (Keq < coefUp*SwapRate  && SwapRate < Keq) ? coefDown*SwapRate: coefDown*Keq;
			sigmaStrike =  Keq < SwapRate + ARM_NumericConstants::ARM_TOLERENCE ?  SwapRate : Keq; 
		}
		/// Built  flat bs model
		double expiry = (swaption->GetExpiryDate().GetJulian() - asOfDate.GetJulian())/K_YEAR_LEN;
		double tenor  = swaption->GetTenor();
		double volatm= oswBSModel->ComputeVol(expiry,tenor,defaultATMK,defaultATMK);

		ARM_VolFlat flatvol(asOfDate,volatm,GetCurrencyUnit());

		flatbsModel.SetVolatility(&flatvol);
		flatbsModel.SetCvxAdjVolatility(&flatvol);
		
		swaption->UpdateStrike(sigmaStrike);
		swaption->SetModel(&flatbsModel);
		prices[0] = swaption->ComputePrice();	
		swaptionList[0]= swaption;	
		vegas[0] = swaption->ComputeSensitivity(K_VEGA)*1.0e-2;
		ARM_StdPortfolioPtr Upportfolio = ARM_StdPortfolioPtr(new ARM_StdPortfolio(swaptionList,weight,prices,vegas));
		
		swaption->UpdateStrike(skewStrike);
		swaption->SetModel(&flatbsModel);
		prices[0] = swaption->ComputePrice();
		swaptionList[0]= swaption;
		vegas[0] = swaption->ComputeSensitivity(K_VEGA)*1.0e-2;
		ARM_StdPortfolioPtr Downportfolio = ARM_StdPortfolioPtr(new ARM_StdPortfolio(swaptionList,weight,prices,vegas));

		/// calib method building
		initVol = volatm*SwapRate/10000.0;
		volCalibMethod.GetCalibParam()->SetValueAtPoint(0,initVol);
		volCalibMethod.SetPortfolio(Upportfolio);

		///Set Only  portfolio
		skewCalibMethod.SetPortfolio(Downportfolio);

		QGmodel.GetModelParams()->SetModelParam(&modelParamVol);
		QGmodel.GetModelParams()->SetModelParam(&modelParamSkew);
		
		/// ATM Calibrate
		skewCalibMethod.Calibrate(&QGmodel);

		double skewValue = ((ARM_CurveModelParam&)(QGmodel.GetModelParams()->GetModelParam(ARM_ModelParamType::Skew))).GetCurve()->GetOrdinates()[0];
		double volValue = ((ARM_CurveModelParam&)(QGmodel.GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility))).GetCurve()->GetOrdinates()[0];


		/// BS Smiled Model
		swaption->SetModel(oswBSModel);
		prices[0] = swaption->ComputePrice();
		swaptionList[0]= swaption;
		vegas[0] = swaption->ComputeSensitivity(K_VEGA)/100.0;

		swaptionpfList[i]= (ARM_Swaption*) swaption->Clone(); 
		pfprices[i] = prices[0];
		pfvegas[i]=vegas[0];

		Downportfolio = ARM_StdPortfolioPtr(new ARM_StdPortfolio(swaptionList,weight,prices,vegas));

		swaption->UpdateStrike(sigmaStrike);
		swaption->SetModel(oswBSModel);
		prices[0] = swaption->ComputePrice();
		swaptionList[0]= swaption;
		vegas[0] = swaption->ComputeSensitivity(K_VEGA)/100.0;
		Upportfolio = ARM_StdPortfolioPtr(new ARM_StdPortfolio(swaptionList,weight,prices,vegas));

		volCalibMethod.GetCalibParam()->SetValueAtPoint(0,volValue);
		volCalibMethod.SetPortfolio(Upportfolio);

		
		static_cast<ARM_CurveModelParam*>(skewCalibMethod.GetCalibParam())->GetUpperBound()[0] = skewValue;
		skewCalibMethod.SetPortfolio(Downportfolio);

		/// Smiled Calibrate
		skewCalibMethod.Calibrate(&QGmodel);

		SkewValues[i] = skewValue = ((ARM_CurveModelParam&)(QGmodel.GetModelParams()->GetModelParam(ARM_ModelParamType::Skew))).GetCurve()->GetOrdinates()[0];
		double skewTime = SkewTimes[i] = ((ARM_CurveModelParam&)(QGmodel.GetModelParams()->GetModelParam(ARM_ModelParamType::Skew))).GetCurve()->GetAbscisses()[0];
		VolValues[i] = volValue = ((ARM_CurveModelParam&)(QGmodel.GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility))).GetCurve()->GetOrdinates()[0];
		double volTime = VolTimes[i]=((ARM_CurveModelParam&)(QGmodel.GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility))).GetCurve()->GetAbscisses()[0];

		weights[i]= Initweight;
		Initweight *= amortweight;
		delete swaption;
		swaption = NULL;
		swaptionpfList[i]->SetModel(NULL);
	}

	itsOWSSkewPF = new ARM_StdPortfolio(swaptionpfList,pfweight,pfprices,pfvegas);
	for(i=0; i<sizepf; i++){
		delete swaptionpfList[i];
		swaptionpfList[i] = NULL;
	}

	///Last step Call Nag to use LSQ Distrubtion
	ARM_GP_Vector precisions(sizepf,SKEW_CALIB_PRECISION);
	ARM_GP_Vector LBoundVector(sizepf,SKEW_LOWER_BOUND);
	ARM_GP_Vector UBoundVector(sizepf,SKEW_UPPER_BOUND);
	QGM_ParameterSet* qgmSet = QGM_CalibrateSkew(&SkewValues,
						&weights,
						&precisions,
						&SkewValues,
						&LBoundVector,
						&UBoundVector,
						0);
    CC_NS(std,auto_ptr)<QGM_ParameterSet> qgmSetPtr(qgmSet);
	ARM_GP_Vector* Skews = qgmSet->get_X();

	//debug
	double objec = qgmSet->get_objective();
	
	///Update Model and CalibMethod
	ARM_GP_Vector LBoundVol(sizepf,SIGMA_MIN_INIT_VALUE);
	ARM_GP_Vector UBoundVol(sizepf,SIGMA_MAX_INIT_VALUE);
	ARM_CurveModelParam  modelParamQ(ARM_ModelParamType::Skew, Skews,&SkewTimes,"SKEW","STEPUPRIGHT");
	ARM_CurveModelParam  modelParamSigma(ARM_ModelParamType::Volatility, &VolValues,&VolTimes,"SIGMA","STEPUPRIGHT",&LBoundVol,&UBoundVol);

	GetPricingModel()->GetModelParams()->SetModelParam(&modelParamQ);
	GetOSWCalibMethod()->SetCalibParam((ARM_CurveModelParam*)modelParamSigma.Clone());
	GetOSWCalibMethod()->SetNbIteration(NB_ITER_MAX);

	if(itsMRSCalibType != ARM_MRSCalibrationType::Unknown)
		GetCalibMethod()->GetNextMethod()->SetNbIteration(NB_ITER_MAX);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: GetOSWCalibMethod
///	Returns: ARM_CalibMethod
///	Action : get calibration method for diagonal swaptions
/////////////////////////////////////////////////////////////////
ARM_CalibMethod* ARM_CRFCalculator::GetOSWCalibMethod() const
{
    if(itsOSWCalibFlag == ARM_SigmaCalibrationType::Unknown)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : No method available because calibration is off");

#ifdef __GP_STRICT_VALIDATION
    if( GetCalibMethod() == ARM_CalibMethodPtr(NULL) ||
        (itsMRSCalibType != ARM_MRSCalibrationType::Unknown && 
		GetCalibMethod()->GetlinkedMethod() == NULL) )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Bootstrap calib method not found");
#endif

    if(itsMRSCalibType != ARM_MRSCalibrationType::Unknown)
        return GetCalibMethod()->GetlinkedMethod();
    else
        return &(*GetCalibMethod());
}
/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: GetOSWBSModel
///	Returns: ARM_BSModel
///	Action : get the BSModel of the Calibration
/////////////////////////////////////////////////////////////////
ARM_BSModel* ARM_CRFCalculator::GetOSWBSModel() const
{
	ARM_BSModel* OSWModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswModelKey]) );
	return OSWModel;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: SetCalibParam
///	Returns: 
///	Action : 
/////////////////////////////////////////////////////////////////
void ARM_CRFCalculator::SetCalibParam(ARM_ModelParam* modelParam)
{
	if(modelParam->GetType() == ARM_ModelParamType::MeanReversion)
	{
		///Re-init the model to have empty params
		ARM_CalibMethodPtr calib2d = GetCalibMethod();
		*(ARM_CurveModelParam*) (GetMktDataManager()->GetData(GetKeys()[MrsKey]))=(modelParam->ToCurveModelParam());
		CreateAndSetModel();
		SetCalibMethod( ARM_CalibMethodPtr( (ARM_CalibMethod*)(*calib2d).Clone() ) );
	}		
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" Model Param Not of a Good Type!");
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: GetSTMCalibMethod
///	Returns: ARM_CalibMethod
///	Action : get calibration method for short term vanillas
/////////////////////////////////////////////////////////////////
ARM_CalibMethod* ARM_CRFCalculator::GetSTMCalibMethod() const
{
    if(itsMRSCalibType == ARM_MRSCalibrationType::Unknown)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : No method available because calibration is off");

#ifdef __GP_STRICT_VALIDATION
    if( GetCalibMethod() == ARM_CalibMethodPtr(NULL))
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Optimisation calib method not found");
#endif

    return &(*GetCalibMethod());
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: GetOSWPortfolio
///	Returns: ARM_Portfolio
///	Action : get the swaption calibration portfolio
/////////////////////////////////////////////////////////////////
const ARM_StdPortfolioPtr ARM_CRFCalculator::GetOSWPortfolio() const
{
    if(itsOSWCalibFlag == ARM_SigmaCalibrationType::Unknown)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : No portfolio available because calibration is off");

#ifdef __GP_STRICT_VALIDATION
    if( GetCalibMethod() == ARM_CalibMethodPtr(NULL) ||
        (itsMRSCalibType != ARM_MRSCalibrationType::Unknown && 
		(GetCalibMethod()->GetlinkedMethod() == NULL ||
            GetCalibMethod()->GetlinkedMethod()->GetPortfolio() == ARM_StdPortfolioPtr(NULL))) )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Bootstrap calib method or portfolio not found");
#endif
    if(itsMRSCalibType != ARM_MRSCalibrationType::Unknown)
        return GetCalibMethod()->GetlinkedMethod()->GetPortfolio();
    else
        return GetCalibMethod()->GetPortfolio();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: SetOSWPortfolio
///	Returns: void
///	Action : affectation of a calibration portfolio made of swaptions
/////////////////////////////////////////////////////////////////
void ARM_CRFCalculator::SetOSWPortfolio(const ARM_StdPortfolio& port)
{
    if(itsOSWCalibFlag == ARM_SigmaCalibrationType::Unknown || !(port.size()))
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : No portfolio settable because calibration is off or Portfolio to set is empty");

     if(port.GetAsset(0)->GetName() != ARM_SWAPTION || !(port.IsSameAssetsName()))
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Diagonal swaption portfolio is not only made of swaptions");

#ifdef __GP_STRICT_VALIDATION
    if( GetCalibMethod() == ARM_CalibMethodPtr(NULL) ||
        (itsMRSCalibType != ARM_MRSCalibrationType::Unknown && 
		GetCalibMethod()->GetlinkedMethod() == NULL) )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Bootstrap calib method not found");
#endif
    
    if(itsMRSCalibType != ARM_MRSCalibrationType::Unknown)
        GetCalibMethod()->GetlinkedMethod()->SetPortfolio( ARM_StdPortfolioPtr((ARM_StdPortfolio*)(const_cast< ARM_StdPortfolio& >(port).Clone())) );
    else
        GetCalibMethod()->SetPortfolio( ARM_StdPortfolioPtr((ARM_StdPortfolio*)(const_cast< ARM_StdPortfolio& >(port).Clone())) );
}
/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: GetSTMPortfolio
///	Returns: ARM_Portfolio
///	Action : get the short term vanilla calibration portfolio
/////////////////////////////////////////////////////////////////
const ARM_StdPortfolioPtr ARM_CRFCalculator::GetSTMPortfolio() const
{

    if(itsMRSCalibType == ARM_MRSCalibrationType::Unknown)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : No portfolio available because calibration is off");

#ifdef __GP_STRICT_VALIDATION
    if( GetCalibMethod() == ARM_CalibMethodPtr(NULL) ||
        GetCalibMethod()->GetPortfolio() == ARM_StdPortfolioPtr(NULL) )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Optimisation calib method or portfolio not found");
#endif

    return GetCalibMethod()->GetPortfolio();
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: SetSTMPortfolio
///	Returns: void
///	Action : affectation of a calibration portfolio made of short
///          term vanilla products
/////////////////////////////////////////////////////////////////
void ARM_CRFCalculator::SetSTMPortfolio(const ARM_StdPortfolio& port)
{
    if(itsMRSCalibType == ARM_MRSCalibrationType::Unknown || !(port.size()))
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" :  No portfolio settable because calibration is off or Portfolio to set is empty");

#ifdef __GP_STRICT_VALIDATION
    if( GetCalibMethod() == ARM_CalibMethodPtr(NULL))
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Optimisation calib method not found");
#endif

    /// Update portfolio
    GetCalibMethod()->SetPortfolio( ARM_StdPortfolioPtr((ARM_StdPortfolio*)(const_cast< ARM_StdPortfolio& >(port).Clone())) );
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: GetCFPortfolio
///	Returns: void
///	Action : affectation of a calibration portfolio made of caplet/floorlet
/////////////////////////////////////////////////////////////////
const ARM_StdPortfolio* ARM_CRFCalculator::GetCFPortfolio() const
{
     if(!itsCapCalibFlag && !itsFloorCalibFlag)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : No portfolio available because Cap & Floor calibrations are off");

    return itsCapFloorPF;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: GetCFPortfolio
///	Returns: void
///	Action : affectation of a calibration portfolio made of caplet/floorlet
/////////////////////////////////////////////////////////////////
const ARM_StdPortfolio* ARM_CRFCalculator::GetSkewPortfolio() const
{
     if(!(itsModelType == ARM_PricingModelType::QGM1F && itsSkewCalFlag))
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : No portfolio available because  the skew calibration is off/model is not QG");

    return itsOWSSkewPF;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: SetCFPortfolio
///	Returns: void
///	Action : affectation of a calibration portfolio made of caplet/floorlet
/////////////////////////////////////////////////////////////////
void ARM_CRFCalculator::SetCFPortfolio(const ARM_StdPortfolio& port)
{
     if(!itsCapCalibFlag && !itsFloorCalibFlag || !(port.size()))
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : No portfolio settable because Cap & Floor calibrations are off");

    /// Test for portfolio consistency
    int i;
	int pfSize = port.size();
    for(i=0;i<pfSize;++i)
    {
        if(port.GetAsset(i)->GetName() != ARM_CAPFLOOR ||
            dynamic_cast< ARM_CapFloor* >(port.GetAsset(i))->GetResetDates()->size() != 1)
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Cap/Floor portfolio is not only made of single caplet/floorlet");
    }

    /// Free memory
    delete itsCapFloorPF;
    DeletePointorVector<ARM_VanillaArg>(itsVanillaArgVect); 

    /// Update portfolio
    itsCapFloorPF = static_cast< ARM_StdPortfolio* >(const_cast< ARM_StdPortfolio& >(port).Clone());
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
    itsVanillaArgVect.resize(pfSize);
    int index;
    for(i=0;i<pfSize;++i)
    {
        index = itsVanillaArgVect[i]->GetIndex();
        delete itsVanillaArgVect[i];
        itsVanillaArgVect[i] = ARM_ConverterFromKernel::ConvertSecuritytoArgObject(port.GetAsset(i),asOfDate);

        /// User must set a CF portfolio consistent with the internal GenSecurity !!
        itsVanillaArgVect[i]->SetIndex(index);
    }
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: SetShortTermCalibFlag
///	Returns: void
///	Action : set the flag to manage the MRS optimisation on
///          short term vanilla. The function maintains the consistency
///          of the CalibMethods object
/////////////////////////////////////////////////////////////////
void ARM_CRFCalculator::SetSTMAndUpdateCalibFlag(ARM_MRSCalibType type)
{
    if(itsMRSCalibType == type)
        return;

    /// Restore calib method and portfolio for diagonal swaptions
    /// (before set the new flag value !)
    ARM_CalibMethod* volCalib = GetOSWCalibMethod();
    ARM_StdPortfolioPtr diagonalSwaptionPF = GetOSWPortfolio();

    itsMRSCalibType = type;

    /// Maintain consistency
    if(itsMRSCalibType != ARM_MRSCalibrationType::Unknown)
    {
        /// Create standard short term vanillas portfolio (MRS calibration)
        ARM_StdPortfolioPtr shortTermVanillaPF(CreateShortTermVanilla(diagonalSwaptionPF));

        /// MRS optimisation with an embedded volatility bootstrapping
        ARM_CalibMethod* mrsCalib = new ARM_CalibMethod(shortTermVanillaPF,ARM_ModelParamVector(),ARM_CalibMethodType::Optimize1D,
                                 MAXITER,ARM_CalibrationTarget::PriceTarget,volCalib);
	    SetCalibMethod(ARM_CalibMethodPtr(mrsCalib));

        /// Compute short term vanilla target prices and MRS initial value(s)
        bool isFreezeWeights=false;
        bool isInitMRS=true;
        ComputeShortTermVanillaPrice(isFreezeWeights,isInitMRS);
    }
    else
    {
        /// STM calib true -> false
        /// Change CalibMethods to be a simple bootstrapping
	    SetCalibMethod(ARM_CalibMethodPtr( (ARM_CalibMethod*)volCalib->Clone() ) );
		/// Update MktDadaManager with Mean Reversion stored in Pricing Model
		ARM_CurveModelParam MRSParam((ARM_CurveModelParam&)(GetPricingModel()->GetModelParams()->GetModelParam(ARM_ModelParamType::MeanReversion)));
		double MRSValue = MRSParam.GetCurve()->GetOrdinates()[0];
		ARM_CurveModelParam* MktmrsParam    = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(GetKeys()[MrsKey]));
		if(MktmrsParam)
			MktmrsParam->SetValueAtPoint(0,MRSValue);
		else
			SetMRS(&MRSParam);
    }
}
/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: CreateAndSetCalibration_Frontier
///	Returns: void
///	Action : create the calibration
/////////////////////////////////////////////////////////////////
void ARM_CRFCalculator::CreateAndSetCalibration_Frontier()
{
	ARM_ModelType   modelType = itsModelType;
	ARM_CurveModelParam  MrsParam = GetMRS();
	if(modelType != ARM_PricingModelType::HWM1F){
		itsModelType = ARM_PricingModelType::HWM1F;
		SetMRS(0.0);
		CreateAndSetModel();
	}
	for ( int k(0); k<itsIsFrontier; ++k){
		Calibrate();

		// Price and frontier
		const ARM_DealDescription& dealDesc = GetGenSecurity()->GetDealDescription();
		size_t colIndex = dealDesc.GetColIndex("Frontier");
		ARM_DealDescriptionPtr subDealDesc = dealDesc.GetSubDescription(1,dealDesc.GetRowsNb()-1,colIndex+1);
		ARM_GenSecurityPtr genSec = ARM_GenSecurityPtr(new ARM_GenSecurity(subDealDesc,GetGenSecurity()->GetPayModelName() ) );

		ARM_GenPricer* genPricer = new ARM_GenPricer( &*genSec,&*GetPricingModel() );
		ARM_AutoCleaner<ARM_GenPricer> HoldGP(genPricer );
		double price = genPricer->Price();

		ARM_VectorPtr frontier = genPricer->GetPricerInfo()->GetContents("Frontier").GetData("Intermediateprices").GetVector();

		// Set Calib at frontier
		ARM_StdPortfolioPtr swaptionPF = GetOSWPortfolio();
		ARM_BSModel* oswBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswModelKey]) );
		for(int i=swaptionPF->size()-1, j = frontier->size()-1; i>=0; --i,--j)
		{
			ARM_Swaption* swopt = static_cast< ARM_Swaption* >(swaptionPF->GetAsset(i));
			double old = swopt->GetStrike();
			swopt->UpdateStrike((*frontier)[i]*100.);
			swopt->SetModel(oswBSModel);
			double price = swopt->ComputePrice();
			double vega = swopt->ComputeSensitivity(K_VEGA)/100.;
			if (vega<VEGA_MIN_TO_REMOVE)
			{
				swopt->UpdateStrike(old);
				price = swopt->ComputePrice();
			}
			swaptionPF->SetPrice(price,i);
		}
	}
	itsModelType=modelType;	
	if(itsModelType != ARM_PricingModelType::HWM1F){
		SetMRS(&MrsParam);
		CreateAndSetModel();
	}

	if(itsModelType == ARM_PricingModelType::QGM1F && itsSkewCalFlag)
        PreCalibrateSkew();
}
/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: CreateAndSetCalibration
///	Returns: void
///	Action : create the calibration
/////////////////////////////////////////////////////////////////
void ARM_CRFCalculator::CreateAndSetCalibration()
{
    /// Create standard diagonal swaption portfolio (sigma calibration)
    ARM_StdPortfolioPtr diagonalSwaptionPF(CreateDiagonalSwaption());

    /// Build a volatility bootstrap calibration on the latter portfolio
    /// The CalibMethod object will be cloned because by default it is not shared
    ARM_CalibMethod volCalib(diagonalSwaptionPF,ARM_ModelParamVector(),ARM_CalibMethodType::Bootstrap1D);

    if(itsMRSCalibType != ARM_MRSCalibrationType::Unknown)
    {
       /// Build a MRS optimisation and embed the volatility bootstrapping
        /// The CalibMethod object will be cloned because by default it is not shared
        ARM_CalibMethod* mrsCalib = new ARM_CalibMethod(ARM_StdPortfolioPtr(NULL),ARM_ModelParamVector(),ARM_CalibMethodType::Optimize1D,
                                 ARM_MAX_ITER,ARM_CalibrationTarget::PriceTarget,&volCalib);

        SetCalibMethod( ARM_CalibMethodPtr(mrsCalib) );

        GetCalibMethod()->SetNextMethod(static_cast< ARM_CalibMethod* >(volCalib.Clone()));
        GetCalibMethod()->GetNextMethod()->SetPortfolio(GetOSWPortfolio());
    }
    else
        SetCalibMethod( ARM_CalibMethodPtr(static_cast< ARM_CalibMethod* >(volCalib.Clone()) ) );

    /// Compute diagonal swaption target prices
    bool isFreezeWeights = false;		// reinit non null weights in the portfolio
    bool isInitSigma	 = true;        // init sigma param for calibration (bounds & init guess)
    bool isUpdateStrike  = true;		// force equivalent strike computation
    if(itsOSWCalibFlag != ARM_SigmaCalibrationType::Unknown)
        ComputeDiagonalSwaptionPrice(isFreezeWeights,isInitSigma,isUpdateStrike);
	
	/// Create standard short term vanillas portfolio (MRS calibration)
	/// and compute short term vanilla target prices
	bool isInitMRS = true; // init MRS param for calibration (bounds & init guess)
	if(itsMRSCalibType == ARM_MRSCalibrationType::stmfirstcolumn)
    {
        ARM_StdPortfolioPtr shortTermVanillaPF(CreateShortTermVanilla(diagonalSwaptionPF));
		GetCalibMethod()->SetPortfolio(shortTermVanillaPF);
		ComputeShortTermVanillaPrice(isFreezeWeights,isInitMRS);
	}
	else if(itsMRSCalibType == ARM_MRSCalibrationType::diagstartfwd)
	{
		ARM_StdPortfolio* shortTermVanillaPF = CreateDiagFWDVanilla(*GetOSWPortfolio());
		GetCalibMethod()->SetPortfolio(ARM_StdPortfolioPtr(shortTermVanillaPF));
		ComputeDiagFWDVanillaPrice(isFreezeWeights,isInitMRS);
	}

	 /// Create the implied caplet & floorlet set
    CreateImpliedCapletFloorlet();

    /// Compute caplet/floorlet target market prices
    ComputeImpliedCapletFloorletPrices(isFreezeWeights);

    if(itsModelType == ARM_PricingModelType::QGM1F && itsSkewCalFlag && !itsIsFrontier)
        PreCalibrateSkew();
   
	if(itsIsFrontier)
		CreateAndSetCalibration_Frontier();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: CreateDiagonalSwaption
///	Returns: a portfolio
///	Action : create the list of diagonal swaptions
/////////////////////////////////////////////////////////////////
ARM_StdPortfolioPtr ARM_CRFCalculator::CreateDiagonalSwaption()
{
    /// Clean standard swap sets
    size_t i;
    for(i=0;i<itsStdSwaps.size();++i)
        delete itsStdSwaps[i].first;
    itsStdSwaps.resize(0);

    for(i=0;i<itsFundStdSwaps.size();++i)
        delete itsFundStdSwaps[i].first;
    itsFundStdSwaps.resize(0);

    for(i=0;i<itsIndexStdSwaps.size();++i)
        delete itsIndexStdSwaps[i].first;
    itsIndexStdSwaps.resize(0);


    ARM_Swap* stdSwap;
	ARM_Swap* fundingSwap;
	ARM_Swap* indexSwap;

    ARM_Swaption* swaption;

    list < ARM_Security* > swaptionList;
    double defaultStrike = -1; // equivalent strike will be calculated later

    const ARM_DealDescription dealDesc = GetGenSecurity()->GetDealDescription();
    size_t nbEvents=dealDesc.GetRowsNb(); // + 1 because of column names

    /// Standard index
    ARM_INDEX_TYPE liborIndex = GetCurrencyUnit()->GetVanillaIndexType();

    /// Get the standard gap of the cpn index
    int stdCpnIndexResetGap = -(GetCurrencyUnit()->GetSpotDays());

    /// Find end date of diagonal swaps
    int fundingEndIdx = GetNextValidDateRow(dealDesc,FundingEndDate,nbEvents-1,-1,0);
	if(fundingEndIdx == 0)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : no funding end date avalaible, please advise!");

    ARM_Date fundSwapEndDate( atof(dealDesc.GetElem(fundingEndIdx,FundingEndDate).c_str()) );
    ARM_Date indexSwapEndDate(fundSwapEndDate);
    ARM_Date stdSwapEndDate(fundSwapEndDate);
    bool isNoStdIndexSwap = itsCpnTiming==K_ARREARS || (itsCpnTiming==K_ADVANCE && itsCpnResetGap != stdCpnIndexResetGap);
    if(isNoStdIndexSwap)
    {
        /// Add one index term to the last index start date to
        /// get the cpn swap end date.
        int lastIndexStartIdx = GetNextValidDateRow(dealDesc,IndexStartDate,nbEvents-1,-1,0);
        char* payCal    = const_cast< char* >(itsCpnPayCal.c_str());
        if(lastIndexStartIdx > 0)
        {
            /// The cpn leg has at least one arrears index i.e it is not always a fixed cpn
            indexSwapEndDate = ARM_Date( atof(dealDesc.GetElem(lastIndexStartIdx,IndexStartDate).c_str()) );
            indexSwapEndDate.AddPeriod(itsCpnIndexTerm);
		    indexSwapEndDate.GoodBusinessDay(K_MOD_FOLLOWING,payCal);
        }
        /// The standard swap end date must cover the last index end date
        /// adding coupons at cpn frequency if necessary
        ARM_Date unadjStdSwapEndDate(stdSwapEndDate);
        while(stdSwapEndDate.GetJulian() < indexSwapEndDate.GetJulian() - ARM_GlobalConstant::ARM_SEVENDAYS_LAG)
        {
            unadjStdSwapEndDate.AddPeriod(itsCpnFreq);
            stdSwapEndDate = unadjStdSwapEndDate;
            stdSwapEndDate.GoodBusinessDay(K_MOD_FOLLOWING,payCal);
        }
    }
	int j=0;
    /// Analyse event dates
    int fundingStartIdx,indexStartIdx;
    for(size_t exerIdx=1;exerIdx<nbEvents;++exerIdx)
    {
        if(atof(dealDesc.GetElem(exerIdx,Fee).c_str()) < NON_CALL_FEE)
        {
            /// Here is a new exercise date, find the corresponding start date for
            /// funding and index coupon legs (pure fixed coupons are skipped)
            fundingStartIdx = GetNextValidDateRow(dealDesc,FundingStartDate,exerIdx,1,nbEvents);
			if(fundingStartIdx == nbEvents)
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : no funding start date avalaible, please advise!");

            ARM_Date fundSwapStartDate( atof(dealDesc.GetElem(fundingStartIdx,FundingStartDate).c_str()) );
            ARM_Date indexSwapStartDate(fundSwapStartDate);

            indexStartIdx = GetNextValidDateRow(dealDesc,IndexStartDate,exerIdx,1,nbEvents);
			/// if it is at least a cpn index
			if(indexStartIdx < nbEvents)                
                indexSwapStartDate = ARM_Date( atof(dealDesc.GetElem(indexStartIdx,IndexStartDate).c_str()) );

            /// Build standard swaps
            fundingSwap=indexSwap=NULL;
            if( fundSwapStartDate.GetJulian() != indexSwapStartDate.GetJulian() ||
                fundSwapEndDate.GetJulian() != indexSwapEndDate.GetJulian()		||
                string(GetCurrencyUnit()->GetCcyName()) != string(GetFundingCcy().GetCcyName()))
            {
                /// In this case build standard swaps corresponding to the funding leg and
                /// to the index coupon leg
                fundingSwap = new ARM_Swap(fundSwapStartDate,fundSwapEndDate,liborIndex,0.0,1.0,K_RCV,
                                            K_DEF_FREQ,K_DEF_FREQ,&GetFundingCcy());
                itsFundStdSwaps.push_back(pair<ARM_Swap*,int>(fundingSwap,fundingStartIdx));

                indexSwap = new ARM_Swap(indexSwapStartDate,indexSwapEndDate,liborIndex,0.0,1.0,K_RCV,
                                            K_DEF_FREQ,K_DEF_FREQ,GetCurrencyUnit());
                itsIndexStdSwaps.push_back(pair<ARM_Swap*,int>(indexSwap,indexStartIdx));

                /// The equivalent standard swap ends at the index standard swap end
                /// (except the case where pure fixed coupons appear at the end of the leg)
                stdSwap = new ARM_Swap(fundSwapStartDate,stdSwapEndDate,liborIndex,0.0,1.0,K_RCV,
										K_DEF_FREQ,K_DEF_FREQ,GetCurrencyUnit());

                /// To save a link to the equivalent standard swap
                itsOtherStdSwapsIdx.push_back(itsFundStdSwaps.size()-1);
            }
            else
            {
                stdSwap  = new ARM_Swap(indexSwapStartDate,indexSwapEndDate,liborIndex,0.0,1.0,K_RCV,
                                        K_DEF_FREQ,K_DEF_FREQ,GetCurrencyUnit());

                itsOtherStdSwapsIdx.push_back(NO_OTHER_STDSWAPS);
            }

            itsStdSwaps.push_back(pair<ARM_Swap*,int>(stdSwap,fundingStartIdx));

            /// Build the market european swaption and set it
            /// (equivalent strike/vol will be computed latter)
            ARM_Date expiryDate((*(stdSwap->GetFloatLeg()->GetResetDates()))[0]);
            swaption = new ARM_Swaption(stdSwap,itsPayRec,K_EUROPEAN,defaultStrike,expiryDate);
            swaptionList.push_back(static_cast< ARM_Security* >(swaption));
			itsKeepOSWIdx.push_back(j++);

        } /// if call date

    } /// for event date

    /// Built OSW portfolio and set default weight
    /// A default price is set to pass the validation test in CalibMethod constructor
    ARM_StdPortfolio* port = new ARM_StdPortfolio(swaptionList);
    for(i=0;i<port->size();++i)
    {
        port->SetWeight(OSW_DEFAULT_WEIGHT,i);
        port->SetPrice((i+1)*OSW_DEFAULT_PRICE,i);
    }

    return ARM_StdPortfolioPtr(port);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: CreateShortTermVanilla
///	Returns: a portfolio
///	Action : create the list of short term vanilla products
/////////////////////////////////////////////////////////////////
ARM_StdPortfolioPtr ARM_CRFCalculator::CreateShortTermVanilla(const ARM_StdPortfolioPtr& diagonalSwaptionPF)
{
    list < ARM_Security* > vanillaList;
    double defaultStrike =-1; // ATM by default and if necessairy the equivalent strike will be calculated later

    const ARM_DealDescription dealDesc = GetGenSecurity()->GetDealDescription();
    size_t nbEvents=dealDesc.GetRowsNb(); // + 1 because of column names


    /// Get the standard fixed leg period and libor index
    int stdFixTermInMonths = 12/GetCurrencyUnit()->GetFixedPayFreq();
    if(stdFixTermInMonths<0)
        stdFixTermInMonths=1;

    string liborTypeName(string(GetCurrencyUnit()->GetCcyName()) == string("EUR") ? "EURIBOR" : "LIBOR");

    string stdLiborTerm=ARM_ArgConvReverse_MatFrequency.GetString(GetCurrencyUnit()->GetLiborTerm());
    string stdLiborTypeName(liborTypeName);
    stdLiborTypeName += stdLiborTerm;
    ARM_INDEX_TYPE stdLiborType = static_cast< ARM_INDEX_TYPE > (ARM_ArgConv_IndexType.GetNumber(stdLiborTypeName));


    /// For each pair of adjacent calibrated diagonal swaption,
    /// associate a short term vanilla product starting at the start of
    /// the longer swaption and ending at the start of the shorter one.
    /// The vanilla product is a swaption if its swap is longer than the
    /// standard fixed leg frequency and a cap otherwise
    size_t i,nbOSW = diagonalSwaptionPF->GetSize();
    ARM_Swaption* swaption;
    ARM_CapFloor* capFloor;
    ARM_Swaption* swaptionStart;
    ARM_Swaption* swaptionEnd;
    ARM_Date startDate,endDate,expiryDate;
    int termInMonths,cfFreq;

	/// To avoid to calibrate same asset 
	int nbSTM = nbOSW -1;
    for(i=0;i <nbSTM;++i)
    {
        swaptionStart=static_cast< ARM_Swaption* >(diagonalSwaptionPF->GetAsset(i));
        startDate=swaptionStart->GetStartDate();
        swaptionEnd=static_cast< ARM_Swaption* >(diagonalSwaptionPF->GetAsset(i+1));
        endDate=swaptionEnd->GetStartDate();
        termInMonths = static_cast<int>(floor(endDate-startDate)/30+0.5);
        if(termInMonths >= stdFixTermInMonths || termInMonths >= 12)
        {
            /// Create a swaption
            ARM_Swap stdSwap(startDate,endDate,stdLiborType,0.0,1.0,K_RCV,
                                    K_DEF_FREQ,K_DEF_FREQ,GetCurrencyUnit());
            expiryDate = (*(stdSwap.GetFloatLeg()->GetResetDates()))[0];
            swaption = new ARM_Swaption(&stdSwap,itsPayRec,K_EUROPEAN,defaultStrike,expiryDate);

            vanillaList.push_back(static_cast< ARM_Security* >(swaption));
        }
		else
        {
            /// Get the nearest standard libor term
            string termStr("1Y");
            cfFreq=K_ANNUAL;
            if(termInMonths<=2)
            {
                termStr="1M";
                cfFreq=K_MONTHLY;
            }
            else if(termInMonths<=4)
            {
                termStr="3M";
                cfFreq=K_QUARTERLY;
            }
            else if(termInMonths<=9)
            {
                termStr="6M";
                cfFreq=K_SEMIANNUAL;
            }


            /// Create a cap/floor with this term
            string cfLiborTypeName(liborTypeName);
            cfLiborTypeName += termStr;
            ARM_INDEX_TYPE cfLiborType = static_cast< ARM_INDEX_TYPE > (ARM_ArgConv_IndexType.GetNumber(cfLiborTypeName));
            capFloor = new ARM_CapFloor(startDate,endDate,(itsPayRec==K_PAY ? K_CAP : K_FLOOR),defaultStrike,
                            cfLiborType,0.0,cfFreq,cfFreq,GetCurrencyUnit());

            vanillaList.push_back(static_cast< ARM_Security* >(capFloor));
        }
    }

    ARM_StdPortfolio* port = new ARM_StdPortfolio(vanillaList);
    for(i=0;i<port->size();++i)
    {
        port->SetWeight(OSW_DEFAULT_WEIGHT,i);
        port->SetPrice((i+1)*OSW_DEFAULT_PRICE,i);
    }

    return ARM_StdPortfolioPtr(port);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: CreateDiagFWDVanilla
///	Returns: a portfolio
///	Action : create the list of Diagonal Swaption forward starting
/////////////////////////////////////////////////////////////////
ARM_StdPortfolio* ARM_CRFCalculator::CreateDiagFWDVanilla(const ARM_StdPortfolio& portfolio)
{
    ARM_StdPortfolio* vanillaPF = (ARM_StdPortfolio*)const_cast<ARM_StdPortfolio&>(portfolio).Clone();
	size_t i, nbVanilla=vanillaPF->GetSize();

    if(itsMRSStrikeType == ARM_MRSStrikeCalibrationType::strikeEquivalent)
    {
	    for(i=0; i< nbVanilla-1; ++i)
	    {
		    double strike = static_cast< ARM_Swaption* >(vanillaPF->GetAsset(i+1))->GetStrike();
		    static_cast< ARM_Swaption* >(vanillaPF->GetAsset(i+1))->UpdateStrike(strike);
	    }
    }
	
	vanillaPF->sortByPrices();
	///To choose only the swaptions most expensive
	ARM_GP_Vector* prices = To_pARM_GP_Vector(vanillaPF->GetMktPrices());
	CC_NS(std,auto_ptr)<ARM_GP_Vector> pricesPtr(prices);
	double sumprices = prices->sum();
	(*prices)/=sumprices;
	int j =0;
	double perc = (*prices)[0];
	while(perc < LEVEL_TO_CHOOSE_PORTFOLIO)
	{
		perc+= (*prices)[j];
		++j;
	}
	ARM_IntVector assetIndex(nbVanilla-j-1);
	for(int k=j+1, p=0; k<nbVanilla; ++k,++p)
		assetIndex[p] = k;
	/// Resize Portfolio
	ARM_StdPortfolio* FinalvanillaPF = vanillaPF-> GetPortfolio(0,j);

	nbVanilla = FinalvanillaPF->size();
	for(i=0; i<nbVanilla; ++i)
    {
		ARM_Swaption* swaption = static_cast< ARM_Swaption* >(FinalvanillaPF->GetAsset(i));

		int FixFrequency = swaption->GetFixedLeg()->GetPaymentFreq();
		int FlaotFrequency = swaption->GetFloatLeg()->GetPaymentFreq();
		if(FixFrequency != FlaotFrequency)
	        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Only swaption with the same frequency is avaliable to create OSWFwd");

		ARM_Vector* PayDates = swaption->GetFloatLeg()->GetPaymentDates();
		size_t floatSize = PayDates->size();

		double firstInterestTerm = swaption->GetFloatLeg()->GetInterestDays()->Elt(0); 

		///test first flow if we have a stub
		int firstindex = 0;
		int firstterm = int(K_YEAR_LEN/firstInterestTerm +0.5);
		ARM_Vector NominalVect(floatSize,100);
		NominalVect[0] = 0.0;
		if(FlaotFrequency != firstterm)
			NominalVect[1] = 0.0;

		ARM_ReferenceValue NominalRef((ARM_Vector*)PayDates->Clone(),(ARM_Vector*)NominalVect.Clone());
		swaption->SetAmount(&NominalRef);
    }

    return FinalvanillaPF;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: GetIndexType
///	Returns: the index type
///	Action : get the ARM index type of the coupon index
/////////////////////////////////////////////////////////////////
ARM_INDEX_TYPE ARM_CRFCalculator::GetIndexType()
{
    string liborTypeName(string(GetCurrencyUnit()->GetCcyName()) == string("EUR") ? "EURIBOR" : "LIBOR");
    string cpnIndexTerm(itsCpnIndexTerm);
    if(cpnIndexTerm=="12M")
        cpnIndexTerm="1Y"; // because we cant convert EURIBOR1Y and not EURIBOR12M
    liborTypeName += cpnIndexTerm;
    return static_cast< ARM_INDEX_TYPE > (ARM_ArgConv_IndexType.GetNumber(liborTypeName));
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: CreateImpliedCapletFloorlet
///	Returns: nothing
///	Action : create the list of the implied caplet & floorlet
/////////////////////////////////////////////////////////////////
void ARM_CRFCalculator::CreateImpliedCapletFloorlet()
{
    const ARM_DealDescription dealDesc = GetGenSecurity()->GetDealDescription();
    size_t nbEvents=dealDesc.GetRowsNb(); // + 1 because of column names

    list< ARM_Security* > capFloorList;
    ARM_CapFloor* capFloor = NULL;

	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();

    /// Get the model name for coupon leg to get a caplet
    /// deal description with the right model name
    string cpnModelName = GetKeys()[itsCpnModelKey];

    size_t nbCF=0;
    double leverage,strike,cpnMin,cpnMax,capStrike,floorStrike;

    /// Build a caplet/floorlet security using the coupon currency
    ARM_INDEX_TYPE liborType = GetIndexType();
    for(size_t eventIdx=1;eventIdx<nbEvents;++eventIdx)
    {
        /// Caplet/Floorlet selection
        if( dealDesc.GetElemFormat(eventIdx,Leverage)		!= ARM_MISSING_TYPE &&
            dealDesc.GetElemFormat(eventIdx,CpnMin)			!= ARM_MISSING_TYPE &&
            dealDesc.GetElemFormat(eventIdx,CpnMax)			!= ARM_MISSING_TYPE &&
            dealDesc.GetElemFormat(eventIdx,IndexStartDate) != ARM_MISSING_TYPE &&
            fabs((leverage = atof(dealDesc.GetElem(eventIdx,Leverage).c_str()))) > K_NEW_DOUBLE_TOL 
			)
        {
            strike      =   atof( dealDesc.GetElem(eventIdx,Strike).c_str() );
            cpnMin      =   atof( dealDesc.GetElem(eventIdx,CpnMin).c_str() );
            cpnMax      =   atof( dealDesc.GetElem(eventIdx,CpnMax).c_str() );

            ARM_Date indexStartDate( atof(dealDesc.GetElem(eventIdx,IndexStartDate).c_str()) );

             /// Compute not adjusted index end date to avoid problem in caplet building
            ARM_Date indexEndDate(indexStartDate);
            indexEndDate.AddPeriod(itsCpnIndexTerm,itsCpnPayCal.c_str());

            /// Build the cap (ARM_Security & VanillaArg versions)
			int CAP = leverage > K_NEW_DOUBLE_TOL ? K_CAP: K_FLOOR;
			int FLOOR = leverage > K_NEW_DOUBLE_TOL ? K_FLOOR: K_CAP;

			double tmpStrike = (strike-cpnMin)/leverage;
            capStrike = tmpStrike >-K_NEW_DOUBLE_TOL ? 100.0*tmpStrike : -K_NEW_DOUBLE_TOL;
            if(capStrike > -K_NEW_DOUBLE_TOL)
            {
                capFloor = new ARM_CapFloor(indexStartDate,indexEndDate,CAP,capStrike,liborType,0.0,K_DEF_FREQ,K_DEF_FREQ,GetCurrencyUnit());
                capFloorList.push_back(static_cast< ARM_Security* >(capFloor));

                /// Convert it to a VanillaArg for further GP pricings...
                /// ...or it could be smart to build and save right now 
                /// the generic security through the CAP keyword !
                itsVanillaArgVect.push_back( ARM_ConverterFromKernel::ConvertSecuritytoArgObject (capFloor,asOfDate.GetJulian(),cpnModelName) );
                itsVanillaArgVect[nbCF]->SetIndex( static_cast<int>(eventIdx) );

                ++nbCF;
            }

            /// If strike is relevant, build the floor (ARM_Security & VanillaArg versions)
           	tmpStrike = (strike-cpnMax)/leverage;
            floorStrike = tmpStrike > -K_NEW_DOUBLE_TOL ? 100.0*tmpStrike : -K_NEW_DOUBLE_TOL;
            if(floorStrike > -K_NEW_DOUBLE_TOL)
            {
                capFloor = new ARM_CapFloor(indexStartDate,indexEndDate,FLOOR,floorStrike,liborType,0.0,K_DEF_FREQ,K_DEF_FREQ,GetCurrencyUnit());
                capFloorList.push_back(static_cast< ARM_Security* >(capFloor));

                /// Convert it to a VanillaArg for further GP pricings...
                /// ...or it could be smart to build and save right now 
                /// the generic security through the CAP keyword !
                itsVanillaArgVect.push_back( ARM_ConverterFromKernel::ConvertSecuritytoArgObject (capFloor,asOfDate.GetJulian(),cpnModelName) );
                itsVanillaArgVect[nbCF]->SetIndex( static_cast<int>(eventIdx) );

                ++nbCF;
            }
			if(capFloor &&capFloor->GetResetDates()->size() != 1)
			{
				CC_Ostringstream os;
				os << ARM_USERNAME << " : Implied Caplet/Floorlet #" << eventIdx << " has more than one flow";
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
			}

        } /// if caplet/floorlet data available

    } /// for event date

    /// Set the new caplet/floorlet portfolio
    delete itsCapFloorPF;
    itsCapFloorPF = new ARM_StdPortfolio(capFloorList);
    for(int i=0;i<itsCapFloorPF->size();++i)
        itsCapFloorPF->SetWeight(CF_DEFAULT_WEIGHT,i);
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: CreateEquivalentStrike
///	Returns: nothing
///	Action : To compute the equivalent strike by taking into account
///          the delta sensitivity of the fixed leg 
/////////////////////////////////////////////////////////////////
double ARM_CRFCalculator::CreateEquivalentStrike(int startRowIdx, int endRowIdx)
{
	ARM_BSModel* cfBSModel = dynamic_cast< ARM_BSModel* >(GetMktDataManager()->GetData(GetKeys()[CfModelKey]) );
    const ARM_DealDescription dealDesc = GetGenSecurity()->GetDealDescription();
    size_t nbEvents=dealDesc.GetRowsNb(); // + 1 because of column names

	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();
    double leverage,cpnMin,cpnMax;
	double fixlegPrice=0.0;
	double SumdeltaitimesLibor = 0.0;
	double sumLibor=0.0;
	double O1Price = 0.0;

	/// Standard index
	ARM_INDEX_TYPE liborIndex = GetCurrencyUnit()->GetVanillaIndexType();
    /// Build a caplet/floorlet security using the coupon currency
    ARM_INDEX_TYPE liborType = GetIndexType();
    for(size_t eventIdx = startRowIdx; eventIdx<=endRowIdx; ++eventIdx)
    {
		double coupon =0.0;
		double deltai =0.0;
		double indexswpPrice = 0.0;
		double strike = atof( dealDesc.GetElem(eventIdx,Strike).c_str() );
        /// Caplet/Floorlet selection
        if( dealDesc.GetElemFormat(eventIdx,Leverage)		!= ARM_MISSING_TYPE &&
            dealDesc.GetElemFormat(eventIdx,CpnMin)			!= ARM_MISSING_TYPE &&
            dealDesc.GetElemFormat(eventIdx,CpnMax)			!= ARM_MISSING_TYPE &&
            dealDesc.GetElemFormat(eventIdx,IndexStartDate) != ARM_MISSING_TYPE &&
            fabs((leverage = atof(dealDesc.GetElem(eventIdx,Leverage).c_str()))) > K_NEW_DOUBLE_TOL 
			)
        {
            cpnMin =   atof( dealDesc.GetElem(eventIdx,CpnMin).c_str() );
            cpnMax =   atof( dealDesc.GetElem(eventIdx,CpnMax).c_str() );

            ARM_Date indexStartDate( atof(dealDesc.GetElem(eventIdx,IndexStartDate).c_str()) );

             /// Compute not adjusted index end date to avoid problem in caplet building
            ARM_Date indexEndDate(indexStartDate);
            indexEndDate.AddPeriod(itsCpnIndexTerm,itsCpnPayCal.c_str());

            /// Build the cap (ARM_Security & VanillaArg versions)
			int CAPFLOOR = leverage > K_NEW_DOUBLE_TOL ? K_FLOOR:K_CAP;
			//int FLOOR = leverage > K_NEW_DOUBLE_TOL ? K_FLOOR: K_CAP;

			double tmpStrike = (strike-cpnMin)/leverage*100;
			tmpStrike = tmpStrike >-K_NEW_DOUBLE_TOL ? tmpStrike : -K_NEW_DOUBLE_TOL;
			ARM_CapFloor cap(indexStartDate,indexEndDate,CAPFLOOR,tmpStrike,liborType,0.0,K_DEF_FREQ,K_DEF_FREQ,GetCurrencyUnit());
			cap.SetModel(cfBSModel);
			double capPrice = cap.ComputePrice()/100.0;
			double capDelta = cap.GetDeltaValues()->Elt(0);

            /// If strike is relevant, build the floor (ARM_Security & VanillaArg versions)
           	tmpStrike = (strike-cpnMax)/leverage *100;
			tmpStrike = tmpStrike > -K_NEW_DOUBLE_TOL ? tmpStrike : -K_NEW_DOUBLE_TOL;
            ARM_CapFloor floor(indexStartDate,indexEndDate,CAPFLOOR,tmpStrike,liborType,0.0,K_DEF_FREQ,K_DEF_FREQ,GetCurrencyUnit());
			floor.SetModel(cfBSModel);
			double floorPrice = floor.ComputePrice()/100.0;
			double floorDelta = floor.GetDeltaValues()->Elt(0);

			double Lvgestrike = 0.0;
			ARM_Swap Indexswaplet(indexStartDate,indexEndDate,liborType,0.0,Lvgestrike,K_PAY,K_DEF_FREQ,K_DEF_FREQ,GetCurrencyUnit());
			Indexswaplet.SetModel(cfBSModel);
			indexswpPrice = Indexswaplet.ComputePrice()/100.0;

			coupon = CAPFLOOR == K_FLOOR ? leverage*(capPrice -floorPrice): leverage*( - capPrice +floorPrice);
			deltai = CAPFLOOR == K_FLOOR ? leverage*(capDelta -floorDelta): leverage*( - capDelta +floorDelta);

        } /// if caplet/floorlet data available
		else
			cpnMin= strike;

		ARM_Date startDate( atof(dealDesc.GetElem(eventIdx,StartDate).c_str()) );

         /// Compute not adjusted index end date to avoid problem in caplet building
        ARM_Date endDate(startDate);
        endDate.AddPeriod(itsCpnIndexTerm,itsCpnPayCal.c_str());

		string cpnIndexTerm(itsCpnIndexTerm);
        int stdFixFreq=K_ANNUAL;
        if(cpnIndexTerm == "1M")
            stdFixFreq=K_MONTHLY;
        else if(cpnIndexTerm =="3M")
            stdFixFreq=K_QUARTERLY;
        else if(cpnIndexTerm == "6M")
            stdFixFreq=K_SEMIANNUAL;

		ARM_FixLeg fixedLeg(startDate,endDate, cpnMin * 100.0,K_RCV, stdFixFreq,itsFixDayCount,K_COMP_PROP,K_ARREARS, K_ADJUSTED,K_SHORTSTART,GetCurrencyUnit());
		fixedLeg.SetModel(cfBSModel);
		coupon+=fixedLeg.ComputePrice()/100.0;

		double margin = atof( dealDesc.GetElem(eventIdx,FundingSpread).c_str() )*100;
		ARM_Swap FundingswapletWithMargin(startDate,endDate,liborIndex,margin,0.0,K_PAY,K_DEF_FREQ,K_DEF_FREQ,GetCurrencyUnit());
		FundingswapletWithMargin.SetModel(cfBSModel);
		double fundWmargin = FundingswapletWithMargin.ComputePrice()/100.0;

		ARM_Swap Fundingswaplet(startDate,endDate,liborIndex,0.0,0.0,K_PAY,K_DEF_FREQ,K_DEF_FREQ,GetCurrencyUnit());
		Fundingswaplet.SetModel(cfBSModel);
		double fundWoutMargin =Fundingswaplet.ComputePrice()/100.0;

		double marginprice = fundWmargin-fundWoutMargin;

		double deltaitimesLibor = deltai*indexswpPrice;
		fixlegPrice += coupon - deltaitimesLibor - marginprice;

		SumdeltaitimesLibor += deltaitimesLibor;
		sumLibor += fundWoutMargin;
		ARM_FixLeg O1(startDate,endDate, 1,K_RCV, stdFixFreq,itsFixDayCount,K_COMP_PROP,K_ARREARS, K_ADJUSTED,K_SHORTSTART,GetCurrencyUnit());
		O1.SetModel(cfBSModel);
		O1Price += O1.ComputePrice()/100;

    } /// for event date

	double Delta = SumdeltaitimesLibor/sumLibor;
	double strikeequi = 1.0/(1.0-Delta)*fixlegPrice /O1Price;

	return strikeequi;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: CreateFloorletAndDelta
///	Returns: nothing
///	Action : create the list of the implied caplet & floorlet
/////////////////////////////////////////////////////////////////
double ARM_CRFCalculator::CreateFloorletAndDelta(int startRowIdx,
     int endRowIdx,const string& evalDateStr)
{
    const ARM_DealDescription dealDesc = GetGenSecurity()->GetDealDescription();
    const string payModelName = GetGenSecurity()->GetPayModelName();

    size_t nbEvents=dealDesc.GetRowsNb(); // + 1 because of column names
	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();

    /// Get the model name for coupon leg to get a caplet
    /// deal description with the right model name
    string cpnModelName = GetKeys()[itsCpnModelKey];
	double leverage;
	
	if(startRowIdx< nbEvents && dealDesc.GetElemFormat(startRowIdx,IndexStartDate) != ARM_MISSING_TYPE)
	{
		ARM_Date indexStartDate( atof(dealDesc.GetElem(startRowIdx,IndexStartDate).c_str()) );
		ARM_Date resetDateDate = indexStartDate.NextBusinessDay(-GetCurrencyUnit()->GetSpotDays(),GetCurrencyUnit()->GetCcyName());
	}

	ARM_Date EndDatefromAs(asOfDate);
	EndDatefromAs.AddPeriod(itsCpnIndexTerm,itsCpnPayCal.c_str());

    /// Build a caplet/floorlet security using the coupon currency
    ARM_INDEX_TYPE liborType = GetIndexType();
	ARM_BSModel* cfBSModel = dynamic_cast< ARM_BSModel* >(GetMktDataManager()->GetData(GetKeys()[CfModelKey]) );
	double sum=0.0;
    double sumtodebug=0.0;
    for(size_t eventIdx=startRowIdx;eventIdx<=endRowIdx;++eventIdx)
    {
        /// Caplet/Floorlet selection
        if( dealDesc.GetElemFormat(eventIdx,Leverage)		!= ARM_MISSING_TYPE &&
            dealDesc.GetElemFormat(eventIdx,CpnMin)			!= ARM_MISSING_TYPE &&
            dealDesc.GetElemFormat(eventIdx,CpnMax)			!= ARM_MISSING_TYPE &&
            dealDesc.GetElemFormat(eventIdx,IndexStartDate) != ARM_MISSING_TYPE &&
            (leverage=atof(dealDesc.GetElem(eventIdx,Leverage).c_str())) > K_NEW_DOUBLE_TOL )
        {
            double strike  =   atof( dealDesc.GetElem(eventIdx,Strike).c_str() );
            double cpnMin  =   atof( dealDesc.GetElem(eventIdx,CpnMin).c_str() );
            double cpnMax  =   atof( dealDesc.GetElem(eventIdx,CpnMax).c_str() );
			double nominal =   atof( dealDesc.GetElem(eventIdx,Nominal).c_str() );

            ARM_Date indexStartDate( atof(dealDesc.GetElem(eventIdx,IndexStartDate).c_str()) );

             /// Compute not adjusted index end date to avoid problem in caplet building
            ARM_Date StartDateFromAs(EndDatefromAs);
            EndDatefromAs.AddPeriod(itsCpnIndexTerm,itsCpnPayCal.c_str());

			ARM_Date indexEndDate(indexStartDate);
            indexEndDate.AddPeriod(itsCpnIndexTerm,itsCpnPayCal.c_str());

            /// Build the cap (ARM_Security & VanillaArg versions)
            double floorStrike = (strike-cpnMin>-K_NEW_DOUBLE_TOL ? 100.0*(strike-cpnMin)/leverage : -K_NEW_DOUBLE_TOL);
            if(floorStrike > -K_NEW_DOUBLE_TOL)
            {
                ARM_CapFloor capFloor(StartDateFromAs,EndDatefromAs,K_CAP,floorStrike,liborType,
                    0.0,K_DEF_FREQ,K_DEF_FREQ,GetCurrencyUnit());

                if(capFloor.GetResetDates()->size() != 1)
                {
	                CC_Ostringstream os;
	                os << ARM_USERNAME << " : Implied Floorlet #" << eventIdx << " has more than one flow";
	                throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
                }
				capFloor.SetModel(cfBSModel);
				double price = capFloor.ComputePrice();
				double vega = capFloor.ComputeSensitivity(K_VEGA)/100;
				vector< ARM_Security* > capFloorList(1,&capFloor);
				vector< double> prices(1,price);
				vector< double> weights(1,1.0);
				vector< double> vegas(1,vega);
				ARM_StdPortfolio* portfolio = new ARM_StdPortfolio(capFloorList,weights,prices,vegas);			


				double vol = capFloor.GetCapletAdjVols()->Elt(0);
				ARM_FlatSurface surface(vol);
				ARM_SurfaceModelParam surfaceModelParam(ARM_ModelParamType::Volatility,&surface,"ImpliedVolatility",0.00,1.0);
				ARM_ModelParamVector calibParams(1,&surfaceModelParam);
				ARM_Normal_ModelParams  nomalParams(calibParams);
				ARM_ZeroCurve* zcCurve = (ARM_ZeroCurve*)cfBSModel->GetZeroCurve()->Clone();
				ARM_Normal_Model model(ARM_ZeroCurvePtr(zcCurve),nomalParams);

                ARM_ModelParamVector sfrmParamVect;
                ARM_CurveModelParam  paramvol(ARM_ModelParamType::Volatility, vol, "Volatililty",0.0, 2.0, TRUE);
                sfrmParamVect.push_back(&paramvol);
                ARM_CurveModelParam  paramMRS(ARM_ModelParamType::MeanReversion, 0.0, "MRS");
                sfrmParamVect.push_back(&paramMRS);
                ARM_CurveModelParam  paramBeta(ARM_ModelParamType::Beta, 1.0, "BETA");
                sfrmParamVect.push_back(&paramBeta);
                
                ARM_INDEX_TYPE liborType = GetIndexType();
                ARM_IRIndex irIndex(liborType,K_DEF_FREQ,K_DEF_FREQ,GetCurrencyUnit());

                ARM_ModelParamsSFRM* sfrmParams = ARM_ModelParamsSFRMFactory.Instance()->CreateModelParamsSFRM(sfrmParamVect,&irIndex,1,K_DIAG);
		        CC_NS(std,auto_ptr)<ARM_ModelParamsSFRM> HoldSFRMParams(sfrmParams);

	            ARM_ZeroCurve* cpnCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));
		        ARM_SFRM sfrmModel( CreateClonedPtr( cpnCurve ),*sfrmParams);

			    ARM_CalibMethod calibmethod(ARM_StdPortfolioPtr(portfolio),calibParams,
				ARM_CalibMethodType::Bootstrap1D, 100, ARM_CalibrationTarget::PriceTarget,NULL,NULL,true);
				calibmethod.Calibrate(&model);

				ARM_CapFloor realCapFloor(indexStartDate,indexEndDate,K_CAP,floorStrike,liborType,
                    0.0,K_DEF_FREQ,K_DEF_FREQ,GetCurrencyUnit());
                
                double capletPrice = GetSubPrice(eventIdx,eventIdx,Caplet,evalDateStr,dealDesc,
                                                        payModelName,ARM_PricingModelPtr((ARM_SFRM*)sfrmModel.Clone()));
                double varLegPrice = GetSubPrice(eventIdx,eventIdx,CpnIndexFlow,evalDateStr,dealDesc,payModelName,itsConvexityModel);
				/// call factory
                price = nominal*ARM_VanillaPricer::Price(&capFloor, &model )/100;
				//price = ARM_VanillaPricer::Price(&realCapFloor, &model )/100;
				sum+= capletPrice;
                sumtodebug+=capletPrice;

                if(capFloor.GetResetDates()->size() != 1)
                {
	                CC_Ostringstream os;
	                os << ARM_USERNAME << " : Implied Floorlet #" << eventIdx << " has more than one flow";
	                throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
                }
            }

        } /// if caplet/floorlet data available
		
    } /// for event date

	return sum;
    //return sumtodebug;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: GetSubPrice
///	Returns: double
///	Action : compute the analytical price of a sub selection of
///          the deal description
/////////////////////////////////////////////////////////////////
double ARM_CRFCalculator::GetSubPrice(int startRowIdx,int endRowIdx,CRFColAlias columnName,
                                      const string& evalDateStr,
                                      const ARM_DealDescription& dealDesc,
                                      const string& payModelName, ARM_PricingModelPtr modelpricing) const
{
    /// Loop for each row between start and end row to set asOfDate because the
    /// GP can't price a GS with all event dates set to a unique date (asOf here)
    double price=0.0;
    ARM_DealDescriptionPtr genDesc;
    ARM_GenSecurity* genLeg;
    ARM_GenPricer* genPricer;
    for(size_t i=startRowIdx;i<=endRowIdx;++i)
    {
	    genDesc=dealDesc.GetSubDescription(i,i,columnName+1);
        genDesc->SetElem(1,0,evalDateStr,ARM_DATE_TYPE);
        genLeg = new ARM_GenSecurity(genDesc,payModelName);
	    genPricer = new ARM_GenPricer( genLeg,&*modelpricing );

        ARM_AutoCleaner< ARM_GenSecurity > HoldGS(genLeg);
        ARM_AutoCleaner< ARM_GenPricer > HoldGP(genPricer);

        /// Compute the analytical leg price (eventDate was changed to asOfDate)
        price += genPricer->Price();
    }
    return price;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: ComputeEquivalentDatas
///	Returns: double
///	Action : compute the equivalent datas of a diagonal swaption
/////////////////////////////////////////////////////////////////
void ARM_CRFCalculator::ComputeEquivalentDatas(
	ARM_Swap* swap,
	ARM_Swaption* swaption,
	ARM_GP_Vector& equivDatas)
{
	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();

    /// Compute the standard swap rate and the equivalent strike
    double swapRate = swap->CptMarketSwapRate();
    double fixedRate = swaption->GetStrike(); // GetDecompStrike() ?

    /// Compute target price & vega
    double price=swaption->ComputePrice();
    double vega=swaption->ComputeSensitivity(K_VEGA);

    ARM_BSModel* oswBSModel = static_cast< ARM_BSModel* >(swaption->GetModel());
    double optMat = (swaption->GetExpiryDate().GetJulian()-asOfDate.GetJulian())/K_YEAR_LEN;
    double swapMat = (swaption->GetEndDate().GetJulian() - swaption->GetStartDate().GetJulian())/K_YEAR_LEN;
    double vol = oswBSModel->ComputeVol(optMat,swapMat,swapRate,fixedRate);

    equivDatas[OSW_SWAP_RATE]=swapRate/100.0;
    equivDatas[OSW_TARGET_VOL]=vol/100.0;
    equivDatas[OSW_TARGET_PRICE]=price;
    equivDatas[OSW_TARGET_VEGA]=vega;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: ComputeEquivalentDatas
///	Returns: double
///	Action : compute the equivalent datas of a diagonal swaption
/////////////////////////////////////////////////////////////////
void ARM_CRFCalculator::ComputeEquivalentDatas(
    const ARM_DealDescription& dealDesc,
    const string& payModelName,
    const pair<ARM_Swap*,int>& stdSwap,
    ARM_Swaption* swaption, 
	const string& evalDateStr, 
	ARM_GP_Vector& equivDatas)
{
    ARM_Swap* swap=stdSwap.first;
    int startIdx = stdSwap.second;

	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();
    size_t nbEvents=dealDesc.GetRowsNb(); /// + 1 because of column names

	double swapRate = swap->CptMarketSwapRate();
	double fixedRate;
	if(itsOSWCalibFlag == ARM_SigmaCalibrationType::strikeATM)
		fixedRate = swapRate;
	
	/// Compute the analytical (eventDate is changed to asOfDate)
	/// floating leg price (funding & cpn index)
	else if(itsOSWCalibFlag == ARM_SigmaCalibrationType::StrikeDeltaApproxi)
	{
		double sstrike = CreateEquivalentStrike(startIdx,nbEvents-1);
		fixedRate = fabs(sstrike);
	}
	else if(itsOSWCalibFlag == ARM_SigmaCalibrationType::strikeEquivalent)
	{
		double varLegPrice = GetSubPrice(startIdx,nbEvents-1,StdVarFlow,evalDateStr,dealDesc,payModelName,itsConvexityModel);
		varLegPrice= fabs(varLegPrice);

		/// Compute the analytical (eventDate is changed to asOfDate)
		/// fixed leg price (strike & funding spread)
		double fixLegPrice = GetSubPrice(startIdx,nbEvents-1,StdFixFlow,evalDateStr,dealDesc,payModelName,itsConvexityModel);

		/// Compute the standard swap rate and the equivalent strike		
		fixedRate = fabs(fixLegPrice)/ varLegPrice * swapRate;
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : no valid type to calibrate sigma");


    /// Take care : use UpdateStrike because SetStrike() only set strike and not
    /// the decomp strike used pricing !!
    swaption->UpdateStrike(fixedRate);

	 /// Get the B&S model from the market data manager
    ARM_BSModel* oswBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswModelKey]) );
	swaption->SetModel(oswBSModel);
	
	/// Compute target price & vega
    double price=swaption->ComputePrice();
    double vega=swaption->ComputeSensitivity(K_VEGA);
    
    double optMat = (swaption->GetExpiryDate().GetJulian()-asOfDate.GetJulian())/K_YEAR_LEN;
    double swapMat = (swaption->GetEndDate().GetJulian() - swaption->GetStartDate().GetJulian())/K_YEAR_LEN;
    double vol = oswBSModel->ComputeVol(optMat,swapMat,swapRate,fixedRate);

    equivDatas[OSW_SWAP_RATE]=swapRate/100.0;
    equivDatas[OSW_TARGET_VOL]=vol/100.0;
    equivDatas[OSW_TARGET_PRICE]=price;
    equivDatas[OSW_TARGET_VEGA]=vega;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: ComputeEquivalentDatas
///	Returns: double
///	Action : compute the equivalent datas of a diagonal swaption
/////////////////////////////////////////////////////////////////
void ARM_CRFCalculator::ComputeEquivalentDatas(
    const ARM_DealDescription& dealDesc,
    const string& payModelName,
    const pair<ARM_Swap*,int>& stdSwap,
    const pair<ARM_Swap*,int>& fundStdSwap,
	const pair<ARM_Swap*,int>& indexStdSwap,
    ARM_Swaption* swaption, 
	const string& evalDateStr, 
	ARM_GP_Vector& equivDatas)
{
	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();
    size_t nbEvents=dealDesc.GetRowsNb(); /// + 1 because of column names

    ARM_BSModel* oswModel = static_cast< ARM_BSModel* >(swaption->GetModel());

    /// Compute the analytical funding leg price (eventDate is changed to asOfDate) and the
    /// standard swap rate corresponding to the cpn index leg
    ARM_Swap* fundingSwap = fundStdSwap.first;
    int fundingStartIdx = fundStdSwap.second;
    double fundVarPrice = GetSubPrice(fundingStartIdx,nbEvents-1,FundingVarFlow,evalDateStr,dealDesc,payModelName,itsConvexityModel);
    double fundSwapRate = fundingSwap->CptMarketSwapRate();

    /// Compute the analytical cpn index leg price (eventDate is changed to asOfDate) and the
    /// standard swap rate corresponding to the cpn index leg
    ARM_Swap* indexSwap = indexStdSwap.first;
    int indexStartIdx = indexStdSwap.second;
	
    double indexPrice = GetSubPrice(indexStartIdx,nbEvents-1,CpnIndexFlow,evalDateStr,dealDesc,payModelName,itsConvexityModel);
	indexPrice = fabs(indexPrice);
    double indexSwapRate = indexSwap->CptMarketSwapRate();

    /// Compute the standard swap rate mixing funding & cpn index legs
    ARM_Swap* swap = stdSwap.first;
    int startIdx = stdSwap.second;
    double swapRate = swap->CptMarketSwapRate();

    /// Compute the analytical fixed leg price (eventDate is changed to asOfDate)
    double fixLegPrice = GetSubPrice(startIdx,nbEvents-1,StdFixFlow,evalDateStr,dealDesc,payModelName,itsConvexityModel);
   
    /// Compute the equivalent strike
    double stdEquivAnnuity = (fabs(fundVarPrice)+indexPrice) / swapRate;
    double fixedRate = fabs(fixLegPrice) / stdEquivAnnuity;

    /// Take care : use UpdateStrike because SetStrike() only set strike and not
    /// the decomp strike used pricing !!
	if(itsOSWCalibFlag == ARM_SigmaCalibrationType::strikeEquivalent)
		swaption->UpdateStrike(fixedRate);

	else if(itsOSWCalibFlag == ARM_SigmaCalibrationType::StrikeDeltaApproxi)
	{
		double sstrike = CreateEquivalentStrike(startIdx,nbEvents-1);
		fixedRate = fabs(sstrike);
		swaption->UpdateStrike(fixedRate);
	}
	else if(itsOSWCalibFlag == ARM_SigmaCalibrationType::strikeATM)
		fixedRate = swapRate;
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : no valid type to calibrate sigma");

    /// Compute the volatility of the standard swap corresponding to the funding leg
    double swapMat = (fundingSwap->GetEndDate().GetJulian() - fundingSwap->GetStartDate().GetJulian())/K_YEAR_LEN;
    double expiryDate = (*(fundingSwap->GetFloatLeg()->GetResetDates()))[0];
    double optMat = (expiryDate-asOfDate.GetJulian())/K_YEAR_LEN;
#ifdef IS_LN_APPROX_K_EQUIV
    double fundVol = oswModel->ComputeVol(optMat,swapMat,fundSwapRate,fixedRate)/100.0;
#else
    double fundVol = oswModel->ComputeVol(optMat,swapMat,fundSwapRate,fixedRate)*fundSwapRate;
#endif

    /// Compute the volatility of the standard swap corresponding to the cpn index leg
    swapMat = (indexSwap->GetEndDate().GetJulian() - indexSwap->GetStartDate().GetJulian())/K_YEAR_LEN;
    expiryDate = (*(indexSwap->GetFloatLeg()->GetResetDates()))[0];
    optMat = (expiryDate-asOfDate.GetJulian())/K_YEAR_LEN;
#ifdef IS_LN_APPROX_K_EQUIV
    double indexVol = oswModel->ComputeVol(optMat,swapMat,indexSwapRate,fixedRate)/100.0;
#else
    double indexVol = oswModel->ComputeVol(optMat,swapMat,indexSwapRate,fixedRate)*indexSwapRate;
#endif

    /// Back to real values
    swapRate /= 100.0;
    fixedRate /= 100.0;
    stdEquivAnnuity *= 100.0;

    /// Compute the equivalent volatility of the swaption (mixing both funding & cpn index leg ones)
    /// through a linear combination with normalized vols
    double fundCoef = fundVarPrice/stdEquivAnnuity;
    double indexCoef = indexPrice/stdEquivAnnuity;
    optMat = (swaption->GetExpiryDate().GetJulian()-asOfDate.GetJulian())/K_YEAR_LEN;

#ifdef IS_LN_APPROX_K_EQUIV
    double stdSwapRateVol = sqrt( log( (fundCoef*fundCoef*exp(fundVol*fundVol*optMat) +
                                        indexCoef*indexCoef*exp(indexVol*indexVol*optMat) +
                                        2*fundCoef*indexCoef*exp(fundVol*indexVol*optMat)) /
                                       (swapRate*swapRate) ) / optMat );
#else
    double stdSwapRateVol = (fundCoef*fundVol + indexCoef*indexVol)/swapRate;
#endif

    /// Compute swaption price & vega
    double price,vega;
    ARM_VolFlat volCurve(const_cast< ARM_Date& >(asOfDate),stdSwapRateVol*100.0); // will be cloned
    ARM_BSModel flatBSModel(oswModel->GetZeroCurve(),&volCurve,K_PRICE);
    swaption->SetModel(&flatBSModel);
    price=swaption->ComputePrice();
    vega=swaption->ComputeSensitivity(K_VEGA);

    /// Reset the local B&S model to disable any implied vol computation
	swaption->SetModelVariable(NULL);			

    equivDatas[OSW_SWAP_RATE]=swapRate;
    equivDatas[OSW_TARGET_VOL]=stdSwapRateVol;
    equivDatas[OSW_TARGET_PRICE]=price;
    equivDatas[OSW_TARGET_VEGA]=vega;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: ReselectSwaption
///	Returns: nothing
///	Action : remove from removable swaptions the most vol sensitive
/////////////////////////////////////////////////////////////////
void ARM_CRFCalculator::ReselectSwaption(ARM_Swaption* swaption,double swapRate,double vol,bool isStd,int idx,ARM_IntVector& removedOSWIdx)
{
	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();
    ARM_BSModel* oswBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswModelKey]) );
	ARM_BSModel flatBSModel;
	if(isStd)
		swaption->SetModel(oswBSModel);
	else
	{
		/// Build an adhoc B&S model
		ARM_VolFlat volCurve(const_cast< ARM_Date& >(asOfDate),vol*100.0);
		flatBSModel = ARM_BSModel(oswBSModel->GetZeroCurve(),&volCurve,K_PRICE);
		swaption->SetModel(&flatBSModel);
	}

    double nominal = swaption->GetAmount()->CptReferenceValue(0.0); // Cst RefValue
	double minVega= VEGA_MIN_TO_REMOVE * (itsMRSCalibType != ARM_MRSCalibrationType::Unknown ? nominal: 1.0);
	double K = swaption->GetStrike();
	swapRate *= 100.0; // %
	double shiftK = 0.05 * (K > swapRate ? -1 : 1); // 5bp
	int i;
	for(i=0;i<100;++i,K+=shiftK)
	{
		swaption->UpdateStrike(K);
		if(swaption->ComputeSensitivity(K_VEGA) > minVega)
			break;
	}

	/// Erase from the list of removable swaption
	ARM_IntVector::iterator iter = removedOSWIdx.find(idx);
	if(iter != removedOSWIdx.end())
		removedOSWIdx.erase(iter);


	if(!isStd)
		/// Reset the adhoc B&S model to disable any implied vol computation
		swaption->SetModelVariable(NULL);			
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: ComputeDiagonalSwaptionPrice
///	Returns: nothing
///	Action : compute market target prices of the swaption portfolio
/////////////////////////////////////////////////////////////////
void ARM_CRFCalculator::ComputeDiagonalSwaptionPrice(bool isFreezeWeights, bool isInitParam, bool isUpdateStrike)
{
    /// Get the B&S model and vol curves from the market data manager
    ARM_BSModel* oswBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswModelKey]) );

    const ARM_DealDescription dealDesc = GetGenSecurity()->GetDealDescription();
    const string payModelName = GetGenSecurity()->GetPayModelName();

	/// Get asOfDate
	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();

    /// Get the standard gap of the cpn index
    char* resetCal = const_cast< char* >(itsCpnResetCal.c_str());
    int stdCpnIndexResetGap = -GetCurrencyUnit()->GetSpotDays();

    if(isUpdateStrike)
    {
        /// Build an SFRM model (instead of a B&S not supported at the
        /// moment by the GP) for convexity ajdustment if necessary
        /// in equivalent strike computation

        ARM_BSModel* cfBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[CfModelKey]) );
	    double defaultValue = 0.0;

        ARM_INDEX_TYPE liborType = GetIndexType();
        ARM_IRIndex irIndex(liborType,K_DEF_FREQ,K_DEF_FREQ,GetCurrencyUnit());

        /// Build a volatility curve with a schedule associated to
        /// the reset date of the coupon index
	    ARM_GP_Vector volTimes;
	    ARM_GP_Vector volValues;

        size_t nbEvents=dealDesc.GetRowsNb(); // + 1 because of column names
        double resetTime,monthTime = K_YEAR_LEN/12.0;
        double liborMat,liborVol;
        double defaultATMK=4.0;
        int fwdRule = GetCurrencyUnit()->GetFwdRule();

        for(size_t eventIdx=1;eventIdx<nbEvents;++eventIdx)
        {
            if(dealDesc.GetElemFormat(eventIdx,IndexStartDate) != ARM_MISSING_TYPE)
            {
                ARM_Date indexStartDate( atof(dealDesc.GetElem(eventIdx,IndexStartDate).c_str()) );
                ARM_Date indexEndDate(indexStartDate);
                indexEndDate.AddPeriod(itsCpnIndexTerm,resetCal);
                indexEndDate.GoodBusinessDay(fwdRule, resetCal);
                ARM_Date indexResetDate(indexStartDate);
                indexResetDate.GapBusinessDay(stdCpnIndexResetGap, resetCal);
                if( (resetTime = indexResetDate.GetJulian() - asOfDate.GetJulian()) > 0.0 )
                {
                    volTimes.push_back(resetTime);
                    liborMat=floor((indexEndDate.GetJulian() - indexStartDate.GetJulian())/monthTime+0.5)/12.0;
                    liborVol = cfBSModel->ComputeVol(resetTime/K_YEAR_LEN,liborMat,defaultATMK,defaultATMK)/100.0;
                    volValues.push_back(liborVol);
                }
            }
        }

        ARM_ModelParamVector sfrmParamVect;
        ARM_CurveModelParam volParam( ARM_CurveModelParam(ARM_ModelParamType::Volatility,&volValues,&volTimes,"SIGMA") );
        sfrmParamVect.push_back(&volParam);
        ARM_CurveModelParam mrsParam( ARM_CurveModelParam(ARM_ModelParamType::MeanReversion,defaultValue,"MRS") );
        sfrmParamVect.push_back(&mrsParam);
        ARM_CurveModelParam shiftParam( ARM_CurveModelParam(ARM_ModelParamType::Shift,defaultValue,"SHIFT") );
        sfrmParamVect.push_back(&shiftParam);

		ARM_ModelParamsSFRM* sfrmParams = 
			ARM_ModelParamsSFRMFactory.Instance()->CreateModelParamsSFRM(sfrmParamVect,&irIndex,1,K_DIAG);
		CC_NS(std,auto_ptr)<ARM_ModelParamsSFRM> HoldSFRMParams(sfrmParams);

	    ARM_ZeroCurve* cpnCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));
		ARM_PricingModel*  sfrmModel = new ARM_SFRM( CreateClonedPtr( cpnCurve ),*sfrmParams);

        if( string(GetCurrencyUnit()->GetCcyName()) != string(GetFundingCcy().GetCcyName()) )
        {
            /// In case of basis, the actual convexity model is the MultiAsset Model
            /// where reference model is the above SFRM model.
            ARM_MultiAssetsModel* newmodel = static_cast< ARM_MultiAssetsModel* >(GetPricingModel()->Clone());

            ///Update the linked ModelMap;
            ARM_ModelNameMap* modelMap = newmodel->GetModelMap();
            ARM_PricingModelPtr oldRefModel = (*modelMap)[myRefModel]->Model();

	        sfrmModel->SetModelName(oldRefModel->GetModelName());
            (*modelMap)[myRefModel]->Model() = ARM_PricingModelPtr(sfrmModel);

            ARM_ForwardMargin* irMargin = static_cast< ARM_ForwardMargin* >(&*((*modelMap)[myIrMarginModel]->Model()));
            irMargin->SetRefModel(sfrmModel);		/// linked model is always the reference model

	        /// compared to update links, we need to know if we need to chang the basis Margin model reference model
            ARM_ForwardMargin* basisMargin = static_cast< ARM_ForwardMargin* >(&*((*modelMap)[myBasisMarginModel]->Model()));
		/*	ARM_ZeroCurve* pCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));
			basisMargin->SetZeroCurve(CreateClonedPtr( pCurve ));*/
            if(basisMargin->GetRefModel() == &*oldRefModel)
                basisMargin->SetRefModel(sfrmModel);

            itsConvexityModel = ARM_PricingModelPtr( newmodel );
        }
        else
            itsConvexityModel = ARM_PricingModelPtr(sfrmModel);
    }
    else
        /// No equivalent strike computation
        itsConvexityModel = GetPricingModel();

    CC_Ostringstream evalDateDesc;
    evalDateDesc << CC_NS(std,fixed) << asOfDate.GetJulian();
    string evalDateStr(evalDateDesc.str());

    ARM_YCModel* ycModel = oswBSModel->GetYCModel();

    ARM_CalibMethod* oswCalibMethod = GetOSWCalibMethod();
    ARM_StdPortfolioPtr swaptionPF = GetOSWPortfolio();
    size_t i,nbOSW=swaptionPF->GetSize();

    /// Test if volatility must be initialise. In this case,
    /// its schedule contents the expiry date of each swaption product
    size_t sigmaIdx,paramSize=oswCalibMethod->GetCalibParams().size();
    for(sigmaIdx=0;sigmaIdx<paramSize;++sigmaIdx)
        if(oswCalibMethod->GetCalibParam(sigmaIdx) &&
           (oswCalibMethod->GetCalibParams())[sigmaIdx]->GetType() == ARM_ModelParamType::Volatility)
            break;

    bool isInitSigma = isInitParam || sigmaIdx >= paramSize || paramSize == 0;

    double price,vega,swapRate,vol,nominal,weight;
    ARM_GP_Vector initTimes(nbOSW);
    ARM_GP_Vector initSigmas(nbOSW);
    ARM_Swaption* swaption;
    ARM_GP_Vector equivDatas(OSW_NB_EQUIVDATA);
    int otherStdSwapsIdx;
    bool isNextCalib = GetCalibMethod()->GetNextMethod() != NULL;
    double initNorVol;
    ARM_IntVector selectedOSWIdx, removedOSWIdx;

	double maxSwapRate,maxVol,maxVega = -ARM_NumericConstants::ARM_BIGGEST_POSITIVE_NUMBER;
	int maxIdx;
	bool maxIsStd,isStd;
    for(int j=0;j<itsKeepOSWIdx.size();++j)
    {
		i = itsKeepOSWIdx[j];
        swaption=static_cast< ARM_Swaption* >(swaptionPF->GetAsset(j));
	    swaption->SetModel(oswBSModel);
        itsStdSwaps[i].first->SetModel(ycModel);
		isStd = true;
        if(!isUpdateStrike)
            /// Just compute price, vega,... but the strike used is already set in the swaption
            ComputeEquivalentDatas(itsStdSwaps[i].first,swaption,equivDatas);
        else if((otherStdSwapsIdx=itsOtherStdSwapsIdx[i]) != NO_OTHER_STDSWAPS)
        {
            /// Compute an equivalent standard swap rate mixing two unsynchronised swaps
            itsFundStdSwaps[otherStdSwapsIdx].first->SetModel(ycModel);
            itsIndexStdSwaps[otherStdSwapsIdx].first->SetModel(ycModel);
            ComputeEquivalentDatas(dealDesc,payModelName,itsStdSwaps[i],itsFundStdSwaps[otherStdSwapsIdx],
                itsIndexStdSwaps[otherStdSwapsIdx],swaption,evalDateStr,equivDatas);
			isStd = false;
        }
        else
            /// Compute an equivalent standard swap rate
            ComputeEquivalentDatas(dealDesc,payModelName,itsStdSwaps[i],swaption,evalDateStr,equivDatas);

        swapRate    = equivDatas[OSW_SWAP_RATE];
        vol         = equivDatas[OSW_TARGET_VOL];
        price       = equivDatas[OSW_TARGET_PRICE];
        vega        = equivDatas[OSW_TARGET_VEGA];
		itsvStrike.push_back(swaption->GetStrike());

		if(vega > maxVega)
		{
			/// Save diagonal swaption with maximum vega
			maxSwapRate	= swapRate;
			maxVol		= vol;
			maxVega		= vega;
			maxIsStd	= isStd;
			maxIdx		= j;
		}

        nominal = swaption->GetAmount()->CptReferenceValue(0.0); // Cst RefValue
        weight=(*(swaptionPF->GetWeights()))[j];
		swaptionPF->SetPrecision(0.001*vega,j);
        swaptionPF->SetPrice(price,j);

		if(isInitSigma)
        {
			double notional = itsMRSCalibType != ARM_MRSCalibrationType::Unknown ? nominal: 1.0;
			bool toRemove = vega < notional*VEGA_MIN_TO_REMOVE;
			if(toRemove)
				removedOSWIdx.push_back(j);
			else if (j == 0)
				selectedOSWIdx.push_back(j);
			else
			{
				bool toSelected = vega < VEGA_MIN_TO_SELECT * nominal && !isFreezeWeights;
				if(toSelected) selectedOSWIdx.push_back(j);
			}
        
            /// Sigma initialisation
            initTimes[j] = swaption->GetExpiryDate().GetJulian()-asOfDate.GetJulian();
            initNorVol = vol*swapRate;
            initSigmas[j] =initNorVol < SIGMA_MIN_INIT_VALUE ? SIGMA_MIN_INIT_VALUE : initNorVol > SIGMA_MAX_INIT_VALUE ? SIGMA_MAX_INIT_VALUE: initNorVol;
        }
    }
	if(isInitSigma)
    {
		/// Check Portfolio
		CC_NS(std,auto_ptr)<ARM_CurveModelParam> modelparam(CheckPortfolioOSW(selectedOSWIdx));

		removedOSWIdx.fill(selectedOSWIdx);
		removedOSWIdx.sort();
		/// Delete unselected swaptions in portfolio and initial curve
		if(removedOSWIdx.size()>0 && removedOSWIdx.size() == nbOSW)
		{
			if(!isUpdateStrike)
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
					" : no swaption with significant vega in calibration portfolio and strikes can't be updated !");

			/// All swaptions were removed because vegas are too low !
			/// => reinsert the maximum vega swaption. Its strike
			/// is shifted so that the swaption is selectable
			swaption=static_cast< ARM_Swaption* >(swaptionPF->GetAsset(maxIdx));
			ReselectSwaption(swaption,swapRate,maxVol,maxIsStd,maxIdx,removedOSWIdx);

		}
		if(removedOSWIdx.size()>0)
		{
			int newNbOSW = nbOSW - removedOSWIdx.size();
			ARM_GP_Vector newInitTimes(newNbOSW),newInitSigmas(newNbOSW);
			ARM_IntVector newKeepOSWIdx(newNbOSW);
			ARM_GP_Vector newStrike(newNbOSW);
			ARM_IntVector::iterator Iter;
			size_t j;
			for(i=0,j=0;i<nbOSW;++i)
			{
				Iter=removedOSWIdx.find(i);
				if(Iter == removedOSWIdx.end())
				{
					newInitTimes[j]    = initTimes[i];
					newInitSigmas[j]   = initSigmas[i];
					newKeepOSWIdx[j]   = itsKeepOSWIdx[i];
					newStrike[j]   = itsvStrike[i];
					++j;
				}
			}
			initTimes     = newInitTimes;
			initSigmas    = newInitSigmas;
			itsKeepOSWIdx = newKeepOSWIdx;
			itsvStrike    = newStrike;
			nbOSW       = newNbOSW;
    	}
		swaptionPF->RemoveAssets(removedOSWIdx);
		if(itsMRSCalibType != ARM_MRSCalibrationType::Unknown &&
			GetCalibMethod()->GetPortfolio() != ARM_StdPortfolioPtr(NULL))
			GetSTMPortfolio()->RemoveAssets(removedOSWIdx);

        /// Replace the sigma param with new initialisations
        ARM_GP_Vector sigmaLowerBound(nbOSW,SIGMA_LOWER_BOUND);
        ARM_GP_Vector sigmaUpperBound(nbOSW,SIGMA_UPPER_BOUND);
        ARM_CurveModelParam* sigma= new ARM_CurveModelParam(ARM_ModelParamType::Volatility,&initSigmas,&initTimes,
            "SIGMA","STEPUPRIGHT",&sigmaLowerBound,&sigmaUpperBound);

		sigma->MergeModelParam(&*modelparam);
        if(sigmaIdx >= paramSize || paramSize == 0)
        {
            oswCalibMethod->GetCalibParams().push_back(sigma);
            if(isNextCalib)
                GetCalibMethod()->GetNextMethod()->GetCalibParams().push_back( static_cast<ARM_CurveModelParam*>(sigma->Clone()) );
        }
        else
        {
            delete oswCalibMethod->GetCalibParam(sigmaIdx);
            (oswCalibMethod->GetCalibParams())[sigmaIdx] = sigma;
            if(isNextCalib)
            {
                delete (GetCalibMethod()->GetNextMethod()->GetCalibParams())[sigmaIdx];
                (GetCalibMethod()->GetNextMethod()->GetCalibParams())[sigmaIdx] = static_cast<ARM_CurveModelParam*>(sigma->Clone());
            }
        }
    }
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: CheckPortfolioOSW
///	Returns: bool
///	Action : to test if any asset in not vega sensitive
/// and use dichtomy method to initialize sigme is necessairy
/////////////////////////////////////////////////////////////////
ARM_CurveModelParam* ARM_CRFCalculator::CheckPortfolioOSW( ARM_IntVector& removedOSWIdx )
{
	size_t sizeRemove = removedOSWIdx.size();
	ARM_StdPortfolioPtr swaptionPF = GetOSWPortfolio();
	ARM_GP_Vector RepairedSigma,times;
	ARM_IntVector removedOSWIdxBis; 
	double time,value =0.0; 
	
	for (int i=0; i< sizeRemove; ++i)
	{
		int j = removedOSWIdx[i];
		ARM_StdPortfolioPtr Pf = ARM_StdPortfolioPtr(swaptionPF->GetPortfolio(j,j+1));		
		ARM_CurveModelParam modelparam(ARM_ModelParamType::Volatility, value,"SIGMA",SIGMA_MIN_INIT_VALUE,SIGMA_MAX_INIT_VALUE );
		ARM_ModelParamVector modelparamVect(1,&modelparam);
		ARM_ModelFitterDes modelfitter(ARM_ModelFitterSolverType::Dichotomy);
		ARM_CalibMethod volFinalCalib(Pf,modelparamVect,ARM_CalibMethodType::Bootstrap1D,&modelfitter,
										ARM_CalibrationTarget::PriceTarget,NULL,NULL,true);

		///Reintialise the local  model and calibrate point by point
		GetPricingModel()->GetModelParams()->SetModelParam(&modelparam);
		volFinalCalib.Calibrate(&*(GetPricingModel()));

		ARM_Bootstrap1D* modelFitter = dynamic_cast< ARM_Bootstrap1D*>(&*(volFinalCalib.GetModelFitter()));
		if(modelFitter->GetSolver()->GetRootState())
		{
			value = modelFitter->GetCalibParam(ARM_ModelParamType::Volatility)->GetValueAtPoint(0);
			time = ((ARM_CurveModelParam*)(modelFitter->GetCalibParam(ARM_ModelParamType::Volatility)))->GetCurve()->GetAbscisses()[0];
			RepairedSigma.push_back(value);
			times.push_back(time);
		}
		else
			removedOSWIdxBis.push_back(j);

	}
	ARM_CurveModelParam modelparam(ARM_ModelParamType::Volatility, 0.0,"SIGMA");
	ARM_ModelParamVector modelparamVect(1,&modelparam);
	
	/// Update the final Model in order to merge it after.
	GetPricingModel()->GetModelParams()->SetModelParam(&modelparam);

	removedOSWIdx.clear();
	removedOSWIdx = removedOSWIdxBis;
	return new ARM_CurveModelParam(ARM_ModelParamType::Volatility,&RepairedSigma,&times,"REPAIREDSIGMA");
}
/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: ComputeShortTermVanillaPrice
///	Returns: nothing
///	Action : compute market target prices of the short term
///          vanilla portfolio
/////////////////////////////////////////////////////////////////
void ARM_CRFCalculator::ComputeShortTermVanillaPrice(bool isFreezeWeights, bool isInitParam,bool isUpdateStrike)
{
    /// Get the B&S models from the market data manager
    ARM_BSModel* oswBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswModelKey]) );
    ARM_BSModel* cfBSModel = dynamic_cast< ARM_BSModel* >(GetMktDataManager()->GetData(GetKeys()[CfModelKey]) );
    ARM_ZeroCurve* cpnCurve     = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));
    ARM_YCModel Ycmodel(cpnCurve);

    ARM_StdPortfolioPtr diagonalSwaptionPF = GetOSWPortfolio();
    size_t i,nbOSW=diagonalSwaptionPF->GetSize();

    ARM_StdPortfolioPtr vanillaPF = GetSTMPortfolio();
    size_t nbVanilla=vanillaPF->GetSize();

    int offset=0;
    if(nbVanilla+1 > nbOSW)
    {
        offset = nbVanilla+1-nbOSW;
    }

    /// Test if MRS must be initialise. In this case and if
    /// MRS is not constant, its schedule contents the expiry
    /// date of each vanilla product
    ARM_CalibMethod* stmCalibMethod = GetSTMCalibMethod();
    size_t mrsIdx,paramSize=stmCalibMethod->GetCalibParams().size();
    for(mrsIdx=0;mrsIdx<paramSize;++mrsIdx)
    {
        if( stmCalibMethod->GetCalibParam(mrsIdx) &&
            (stmCalibMethod->GetCalibParams())[mrsIdx]->GetType() == ARM_ModelParamType::MeanReversion )
            break;
    }

    bool isInitMRS = isInitParam || mrsIdx >= paramSize || paramSize == 0;

    ARM_CurveModelParam* mrsParam = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(GetKeys()[MrsKey]));
    size_t nbMRS = mrsParam->size();
	if(nbMRS > 1) stmCalibMethod->SetMethodType(ARM_CalibMethodType::Optimize);
    if(isInitMRS && nbMRS > 1)
        nbMRS=nbVanilla;
    ARM_GP_Vector initTimes(nbMRS,0.0);
    ARM_GP_Vector initValues(nbMRS,MRS_DEFAULT_VALUE);

    ARM_Swaption* swaption;
    ARM_CapFloor* capFloor;
    ARM_Security* vanilla;
    double precision,minPrec,price,vega,nominal,weight;

    /// Restore equivalent strikes of diagonal swaption portfolio then compute
    /// target prices and vegas for calibration precision
	/// Get asOfDate
    double defaultATMK = 2.50;
	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();
    bool isNotSelected;
    int idx;
    for(i=0;i<nbVanilla;++i)
    {
		idx = i<=offset ? 0 : i-offset;

        vanilla = vanillaPF->GetAsset(i);
        if(isInitMRS && nbMRS > 1)
            initTimes[i] = vanilla->GetExpiryDate().GetJulian()-asOfDate.GetJulian();

        swaption = dynamic_cast< ARM_Swaption* >(vanilla);
        if(swaption)
        {
			/// Warning, don't neither delete nor move it
			/// else ARM service for Summit will crash
			swaption->SetModel(oswBSModel);
			if(isUpdateStrike)
			{

				double strike = static_cast< ARM_Swaption* >(diagonalSwaptionPF->GetAsset(idx))->GetStrike();
                if(itsMRSStrikeType == ARM_MRSStrikeCalibrationType::strikeATM)
					strike = -1.0;

                else if(itsMRSStrikeType == ARM_MRSStrikeCalibrationType::strikeKeepStdMoyenessCst)
                {
                    static_cast< ARM_Swaption* >(diagonalSwaptionPF->GetAsset(idx))->SetModel(oswBSModel);
                    double fwdDiag = static_cast< ARM_Swaption* >(diagonalSwaptionPF->GetAsset(idx))->PriceToRate(const_cast<ARM_Date&> (asOfDate),0.0);
                    double expiry = (static_cast< ARM_Swaption* >(diagonalSwaptionPF->GetAsset(idx))->GetExpiryDate().GetJulian() - asOfDate.GetJulian())/K_YEAR_LEN;
                    double tenor = static_cast< ARM_Swaption* >(diagonalSwaptionPF->GetAsset(idx))->GetTenor();
		            double oswvolatm= oswBSModel->ComputeVol(expiry,tenor,defaultATMK,defaultATMK);
                    
                    swaption->SetModel(oswBSModel);
                    double fwd = swaption->PriceToRate(const_cast<ARM_Date&> (asOfDate),0.0);
                    tenor  = swaption->GetTenor();
                    double stmvolatm= oswBSModel->ComputeVol(expiry,tenor,defaultATMK,defaultATMK);

                    double tmpstrike = (strike - fwdDiag)/fwdDiag/oswvolatm * fwd*stmvolatm + fwd;
                    strike = tmpstrike < K_NEW_DOUBLE_TOL ? K_NEW_DOUBLE_TOL : tmpstrike;
                }
				swaption->UpdateStrike(strike);
			}
        }
        else
        {
            capFloor = dynamic_cast< ARM_CapFloor* >(vanilla);
            if(capFloor)
            {
                double strike = static_cast< ARM_Swaption* >(diagonalSwaptionPF->GetAsset(idx))->GetStrike();

                capFloor->SetStrike(strike);
                capFloor->SetModel(cfBSModel);
            }

#ifdef __GP_STRICT_VALIDATION
            else
                ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : short term vanilla portfolio has other products that swaption & cap/floor");
#endif
        }

        /// Compute price & vega
        price=vanilla->ComputePrice();
        vega=vanilla->ComputeSensitivity(K_VEGA);

        nominal = vanilla->GetAmount()->CptReferenceValue(0.0); // Cst RefValue
        weight=(*(vanillaPF->GetWeights()))[i];

        isNotSelected = vega < VEGA_MIN_TO_SELECT * nominal && !isFreezeWeights;

        vanillaPF->SetWeight(isNotSelected ? 0.0 : weight,i);

        /// Set precision to max(1/4 vega,0.5bp)
		if(string(GetCurrencyUnit()->GetCcyName()) == string("JPY"))
			vanillaPF->SetPrecision(0.00001,i);
		else{
			precision = 0.25*vega;
			minPrec = 0.00005*nominal;
			if(precision < minPrec)
				precision = minPrec;
			vanillaPF->SetPrecision(precision,i);
		}
        vanillaPF->SetPrice(price,i);
    }

	if(vanillaPF->size() ==0 )
	{

		itsMRSCalibType = ARM_MRSCalibrationType::Unknown;
		ARM_CalibMethod* sigCalibMethod = static_cast<ARM_CalibMethod*>(GetCalibMethod()->GetNextMethod()->Clone());
		SetCalibMethod( ARM_CalibMethodPtr(sigCalibMethod) );
		return;
	
	}
    if(isInitMRS)
    {
        ARM_Curve initMRS(initTimes,initValues,new ARM_StepUpRightOpenCstExtrapolDble);

        if(nbVanilla == 0)
        {
            /// No more than one exercise in the CRF => MRS calibration will do
            /// nothing and will return the MRS initial value => set it to the input MRS
            initMRS = *(mrsParam->GetCurve());
        }
		
        ARM_GP_Vector mrsLowerBound;
        ARM_GP_Vector mrsUpperBound;

		mrsLowerBound = ARM_GP_Vector(nbMRS,MRS_LOWER_BOUND);
		mrsUpperBound = ARM_GP_Vector(nbMRS,MRS_UPPER_BOUND);

        ARM_CurveModelParam* mrs = new ARM_CurveModelParam( ARM_ModelParamType::MeanReversion,
			&initMRS.GetOrdinates(), &initMRS.GetAbscisses(), "MRS", "STEPUPRIGHT", &mrsLowerBound, &mrsUpperBound);

        if(mrsIdx >= paramSize || paramSize == 0)
            stmCalibMethod->GetCalibParams().push_back(mrs);
        else
        {
            delete stmCalibMethod->GetCalibParam(mrsIdx);
            (stmCalibMethod->GetCalibParams())[mrsIdx] = mrs;
        }
    }
}
/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: ComputeDiagFWDVanillaPrice
///	Returns: nothing
///	Action : compute market target prices of the short term
///          vanilla portfolio
/////////////////////////////////////////////////////////////////
void ARM_CRFCalculator::ComputeDiagFWDVanillaPrice(bool isFreezeWeights, bool isInitParam)
{
    /// Get the B&S models from the market data manager
	ARM_BSModel* oswBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswModelKey]) );

    ARM_StdPortfolioPtr vanillaPF = GetSTMPortfolio();
	size_t k, nbVanilla=vanillaPF->GetSize();

    /// Test if MRS must be initialise. In this case and if
    /// MRS is not constant, its schedule contents the expiry
    /// date of each vanilla product
    ARM_CalibMethod* stmCalibMethod = GetSTMCalibMethod();

    

    ARM_Swaption* swaption;
    double price;
	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();
	vector< string > keys(2,"YC_JPY");
	keys[1]="YC_OWSMOD";

	ARM_MarketData_ManagerRep mktDataManager(asOfDate);
	ARM_ZeroCurve* cpnCurve     = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));
	mktDataManager.RegisterData("YC_JPY",cpnCurve);
	mktDataManager.RegisterData("YC_OWSMOD",oswBSModel);
	ARM_MarketIRModel mktIrModel(mktDataManager, keys);

    /// Restore equivalent strikes of diagonal swaption portfolio then compute
    /// target prices and vegas for calibration precision
	/// Get asOfDate
	ARM_YCModel Ycmodel(cpnCurve);
    for(k=0; k<nbVanilla; ++k)
    {
		swaption = dynamic_cast< ARM_Swaption* >(vanillaPF->GetAsset(k));
		swaption->SetModel(&Ycmodel);
		if(itsMRSStrikeType == ARM_MRSStrikeCalibrationType::strikeATM)
		{
			double fwd = swaption->PriceToRate(const_cast<ARM_Date&> (asOfDate),0.0);
			swaption->UpdateStrike(fwd);
			swaption->SetModel(oswBSModel);
		}

        /// Compute price & vega
        //double vega=swaption->ComputeSensitivity(K_VEGA);
		price =	ARM_VanillaPricer::Price(swaption,& mktIrModel );

        /// Set precision to max(1/4 vega,0.5bp)
       // double precision = 0.25*vega;
        vanillaPF->SetPrecision(1.0e-4,k);

        vanillaPF->SetPrice(price,k);
		swaption->SetModel(NULL);
   }

	ARM_CurveModelParam* mrsParam = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(GetKeys()[MrsKey]));

	ARM_GP_Vector lowerBound(1,MRS_LOWER_BOUND);
	ARM_GP_Vector  upperBound(1,MRS_UPPER_BOUND);
	mrsParam->SetLowerBound((ARM_GP_Vector*)lowerBound.Clone());
	mrsParam->SetUpperBound((ARM_GP_Vector*)upperBound.Clone());

	stmCalibMethod->GetCalibParams().push_back(mrsParam);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: ComputeImpliedCapletFloorletPrices
///	Returns: nothing
///	Action : compute market target prices of the C/F portfolio
/////////////////////////////////////////////////////////////////
void ARM_CRFCalculator::ComputeImpliedCapletFloorletPrices(bool isFreezeWeights)
{
    /// Get the B&S model and vol curves from the market data manager
    ARM_BSModel* cfBSModel = dynamic_cast< ARM_BSModel* >(GetMktDataManager()->GetData(GetKeys()[CfModelKey]) );

#ifdef __GP_STRICT_VALIDATION
    if(!itsCapFloorPF ||
        itsCapFloorPF->GetSize() != itsVanillaArgVect.size() )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : can't compute implied cap/floor target price because of inconsistency in the portfolio");
#endif

    /// Compute target market prices with an ARM pricing
    /// Save data in the ARM portfoltio (for market target
    /// price recalculation and exportation) and in the VanillaArg
    /// struct for further strike adjustment (with GP pricing !)
    double price,kappa,nominal,weight;
    ARM_CapFloor* capFloor;
    size_t nbCF=itsCapFloorPF->GetSize();
    bool isNotSelected;
    for(size_t i=0;i<nbCF;++i)
    {
        capFloor=static_cast< ARM_CapFloor* >(itsCapFloorPF->GetAsset(i));
	    capFloor->SetModel(cfBSModel);

        kappa=capFloor->ComputeSensitivity(K_KAPPA);
        kappa=fabs(kappa);
        price=capFloor->ComputePrice();

        nominal = capFloor->GetAmount()->CptReferenceValue(0.0); // Cst RefValue
        weight=(*(itsCapFloorPF->GetWeights()))[i];

        isNotSelected = kappa < KAPPA_MIN_TO_SELECT * nominal && !isFreezeWeights;

        itsCapFloorPF->SetWeight(isNotSelected ? 0.0 : weight,i);
        itsCapFloorPF->SetPrecision(0.001*kappa,i);
        itsCapFloorPF->SetPrice(price,i);

        itsVanillaArgVect[i]->SetMktPrice( price );
    }
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: Calibrate (CalibrateStrikeSpread)
///	Returns: nothing
///	Action : perform the calibration of strike adjustments so that
///          model caplet & floorlet prices fit market ones
/////////////////////////////////////////////////////////////////
void ARM_CRFCalculator::Calibrate()
{
    ARM_MultiAssetsModel* model= dynamic_cast< ARM_MultiAssetsModel* >(&*GetPricingModel());

    if(itsOSWCalibFlag != ARM_SigmaCalibrationType::Unknown)
        GetCalibMethod()->Calibrate(&*(*model->GetModelMap())[GetKeys()[YcKey]]->Model());

    size_t nbCF=itsCapFloorPF->GetSize();

#ifdef __GP_STRICT_VALIDATION
    if(nbCF != itsVanillaArgVect.size())
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : insconsistency in auto-calibrated cap/floor portfolio");
#endif

    CRFColAlias spreadCol;
    string zeroValue("0.0");

    /// Make a copy of the initial description to update strike spreads
    ARM_DealDescriptionPtr dealDesc( (ARM_DealDescription*) const_cast< ARM_DealDescription& >(GetGenSecurity()->GetDealDescription()).Clone() );

    if(itsCapCalibFlag || itsFloorCalibFlag)
    {
        ARM_VanillaCapDigitalArg* capFloorArg;
        const double guessRatioMin=0.75;
        const double guessRatioMax=1.125;

        CRFStrikeAdjFct strikeAdjFct(0,this);

        /// Create a local solver
		const double target = 0.0;
		const double ftol = 1.0e-8;
		const double xtol = 1.0e-8;
		const size_t max_iter = 50;
        T_BrentSolver<CRFStrikeAdjFct> Solver(strikeAdjFct,target,ftol,xtol,max_iter);

        double strike,spread;
        for(size_t i=0;i<nbCF;++i)
        {
            if( ( (itsCapCalibFlag && itsVanillaArgVect[i]->GetCallPut() == K_CAP) ||
                  (itsFloorCalibFlag && itsVanillaArgVect[i]->GetCallPut() == K_FLOOR) ) &&
                (*(itsCapFloorPF->GetWeights()))[i] > K_NEW_DOUBLE_TOL)
            {
                strikeAdjFct.SetProdIdx(i);


                /// Convert to CapArg (because of itsStrikes) then update solver bounds
                capFloorArg = dynamic_cast< ARM_VanillaCapDigitalArg* >(itsVanillaArgVect[i]);

#ifdef __GP_STRICT_VALIDATION
                if(!capFloorArg)
                    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : caplet/floorlet expected for strike spread calibration");
#endif
                strike = (*capFloorArg->GetStrikes())[0];

                /// If caplet is to be calibrated test if it is worth something
                double price = strikeAdjFct(strike) + capFloorArg->GetMktPrice();
                if(price < K_NEW_DOUBLE_TOL && fabs(strike ) > K_NEW_DOUBLE_TOL)
                    itsCapCalibFlag && itsVanillaArgVect[i]->GetCallPut() == K_CAP ? Solver.setInitialGuess(0,0.5*strike,strike):
																	                 Solver.setInitialGuess(0,strike,2*strike);
                else 
					fabs(strike)<= K_NEW_DOUBLE_TOL ? Solver.setInitialGuess(0,-0.1,0.1): Solver.setInitialGuess(0,strike*guessRatioMin,strike*guessRatioMax);


                /// Compute the strike spread
                spread = Solver.Solve() - strike;

                /// Restore initial strike (changed by the solver)
				(*(capFloorArg->GetStrikes()))[0]=strike;


                /// Set the strike spread in the deal description
                CC_Ostringstream adjSpreadDesc;
                adjSpreadDesc << CC_NS(std,fixed) << CC_NS(std,setprecision)(8) << spread;
                spreadCol = (itsVanillaArgVect[i]->GetCallPut() == K_CAP ? AdjCapSpread : AdjFloorSpread);
                dealDesc->SetElem(itsVanillaArgVect[i]->GetIndex(),spreadCol,adjSpreadDesc.str(),ARM_DOUBLE);
            }
            else
            {
                spreadCol = (itsVanillaArgVect[i]->GetCallPut() == K_CAP ? AdjCapSpread : AdjFloorSpread);
                dealDesc->SetElem(itsVanillaArgVect[i]->GetIndex(),spreadCol,zeroValue,ARM_DOUBLE);
            }
        }/// for nb CF
    }
    else
    {
        /// Force Strike Spreads to zero
        for(size_t i=0;i<nbCF;++i)
        {
            spreadCol = (itsVanillaArgVect[i]->GetCallPut() == K_CAP ? AdjCapSpread : AdjFloorSpread);
            dealDesc->SetElem(itsVanillaArgVect[i]->GetIndex(),spreadCol,zeroValue,ARM_DOUBLE);
        }
    }

    /// Rebuild the GenSec to update nodes of the syntaxic tree
	SetGenSecurity( ARM_GenSecurityPtr( new ARM_GenSecurity( dealDesc, GetGenSecurity()->GetPayModelName() ) ) );

}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: GetSubGenSecurity
///	Returns: ARM_GenSecurity*
///	Action : ability to get a sub security
/////////////////////////////////////////////////////////////////

ARM_GenSecurity* ARM_CRFCalculator::GetSubGenSecurity(int prodIdx) const 
{
	/// Build the sub-generic security
	size_t startRowIdx = 1; /// start with the first actual row (=> skip column names)
	size_t endRowIdx = GetGenSecurity()->GetDealDescription().GetRowsNb() - 1;
	ARM_DealDescriptionPtr subDealDesc = GetGenSecurity()->GetDealDescription().GetSubDescription(startRowIdx,endRowIdx,prodIdx+1);
	ARM_GenSecurity* genSec = new ARM_GenSecurity(subDealDesc,GetGenSecurity()->GetPayModelName());
	return genSec;
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: Price
///	Returns: a double
///	Action : price the CRF deal. An auto-calibration is done before
///          by bootstrapping the volatility on diagonal swaptions
///          and computing strike spreads of implied caplets/floorlets
/////////////////////////////////////////////////////////////////
double ARM_CRFCalculator::Price()
{
    /// Caplet/floorlet strike spread adjustment
    CalibrateAndTimeIt();

    /// Price the implicit product according to internal flag
    ARM_GenSecurityPtr genSec = GetGenSecurity();

    if(itsProductToPrice != CRF_PRICE)
    {
        /// Select the right column that describes the product
        CRFColAlias prodIdx;

        if(itsProductToPrice == CAP_PRICE)
            prodIdx = Caplet;
        else if(itsProductToPrice == FLOOR_PRICE)
            prodIdx = Floorlet;
        else if(itsProductToPrice == FUNDING_PRICE)
            prodIdx = Funding;
        else if(itsProductToPrice == STDLEG_PRICE)
            prodIdx = StdFlow;
        else if(itsProductToPrice == RFLEG_PRICE)
            prodIdx = RFFlow;
        else if(itsProductToPrice == STDSWAP_PRICE)
            prodIdx = StdSwaplet;
        else if(itsProductToPrice == RFSWAP_PRICE)
            prodIdx = RFSwaplet;
        else if(itsProductToPrice == BERMUDA_PRICE)
            prodIdx = StdBermudaPrice;

        /// Build the sub-generic security
        size_t startRowIdx = 1; /// start with the first actual row (=> skip column names)
        size_t endRowIdx = genSec->GetDealDescription().GetRowsNb() - 1;
        ARM_DealDescriptionPtr subDealDesc = genSec->GetDealDescription().GetSubDescription(startRowIdx,endRowIdx,prodIdx+1);
        genSec = ARM_GenSecurityPtr(new ARM_GenSecurity(subDealDesc,genSec->GetPayModelName()));
    }

    ARM_GenPricer* genPricer = new ARM_GenPricer( &*genSec,&*GetPricingModel() );
    ARM_AutoCleaner<ARM_GenPricer> HoldGP(genPricer );
    double price  = genPricer->Price();
	price *= GetPorS();
    return price;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: ComputePricingData
///	Returns: a double
///	Action : get pricing data of the calculaltor
/////////////////////////////////////////////////////////////////

void ARM_CRFCalculator::ComputePricingData() const
{}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: UpdateModel
///	Returns: void
///	Action : update the model datas w.r.t. MarketDataManager
/////////////////////////////////////////////////////////////////
void ARM_CRFCalculator::UpdateModel()
{
	ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();
	/// Get constant Mean Reversion from Curve Model Param.
	ARM_CurveModelParam* mrsParam = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(GetKeys()[MrsKey]));
	//double  flowStartDate=(*(datesStructure.GetDateStrip(EXER_SCHED)->GetFlowStartDates()))[itsFirstCallIdx];
	double tenorLag = itsStartDate > asOfDate ? itsEndDate - itsStartDate : itsEndDate - asOfDate;
	double mrsValue = mrsParam->GetValue(tenorLag);

	ARM_CurveModelParam CstmrsParam( ARM_ModelParamType::MeanReversion, mrsValue,"MRS");

	SetMRS(&CstmrsParam);

	ARM_ZeroCurve* cpnCurve     = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));
	ARM_ZeroCurve* fundingCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[FundingKey]));
	ARM_ZeroCurve* basisCurve   = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[BasisKey]));

    /// Check for model consistency
    ARM_MultiAssetsModel* model= dynamic_cast< ARM_MultiAssetsModel* >(&*GetPricingModel());
    if( !model)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Hybrid BFIR model is not a good type for updating");

    ARM_ModelNameMap* modelMap = model->GetModelMap();
	/// Allow update only for default model (H&W1F/QGM1F)
    ARM_HullWhite1F* RefModel = dynamic_cast< ARM_HullWhite1F* >(&*(*modelMap)[GetKeys()[YcKey]]->Model());
    if( !RefModel )
	{
		 ARM_QGM1F* RefModel = dynamic_cast< ARM_QGM1F* >(&*(*modelMap)[GetKeys()[YcKey]]->Model());
		 if( !RefModel )
			 ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+": Only HWM1F/QGM1F model are allowed for the stochastic Model");
	}
   
    /// To update the MRS param
    (*modelMap)[GetKeys()[YcKey]]->Model()->GetModelParams()->SetModelParam(&CstmrsParam);

     model->SetZeroCurve(CreateClonedPtr(cpnCurve));
     (*modelMap)[myRefModel]->Model()->SetZeroCurve(CreateClonedPtr(cpnCurve ));

    /// Update yield curves
    if(IsBasis())
    {
	    ARM_Forex* forex   = dynamic_cast<ARM_Forex*>(GetMktDataManager()->GetData(GetKeys()[ForexKey]));
        (*modelMap)[myIrMarginModel]->Model()->SetZeroCurve(CreateClonedPtr(fundingCurve));
        (*modelMap)[myBasisMarginModel]->Model()->SetZeroCurve(CreateClonedPtr(basisCurve));

        static_cast< ARM_ForwardForex&>( *((*modelMap)[myForexModel]->Model()) ).SetForex(*forex);

        /// Set forex domestic curve
        static_cast< ARM_ForwardForex&>( *((*modelMap)[myForexModel]->Model()) ).SetCurves(
            (*modelMap)[(*modelMap)[myForexModel]->OtherModelRefNb()[0]]->Model()->GetZeroCurve(),
            (*modelMap)[(*modelMap)[myForexModel]->OtherModelRefNb()[1]]->Model()->GetZeroCurve() );
    }
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: UpdateCalibration
///	Returns: void
///	Action : update the calibration datas w.r.t. MarketDataManager
///          The context is an hedge ratios computation
/////////////////////////////////////////////////////////////////
void ARM_CRFCalculator::UpdateCalibration(bool isUpdateStrike)
{
    if(itsMRSCalibType != ARM_MRSCalibrationType::Unknown)
    {
        /// Share the OSW portfolio between the linked calib method, the bootstrap on vol,
        /// and the next calib method, the final bootstrap on vol
        /// (as it is done at the calculator creation)
        GetCalibMethod()->GetNextMethod()->SetPortfolio(GetCalibMethod()->GetlinkedMethod()->GetPortfolio());
    }

    /// Get the current market model for swaption and its associated YC model
    ARM_BSModel* oswBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswModelKey]) );

    /// Update diagonal swaption prices (and equivalent strike & vol if necessary)
    bool isFreezeWeights = true; // keep the portfolio size to avoid hedge jumps
    bool isInitSigma = false;       // keep current sigma param to init calibration
    if(itsOSWCalibFlag != ARM_SigmaCalibrationType::Unknown)
        ComputeDiagonalSwaptionPrice(isFreezeWeights,isInitSigma,isUpdateStrike);

	/// Update short term vanilla prices
        bool isInitMRS = false; // keep current MRS param to init calibration
	if(itsMRSCalibType == ARM_MRSCalibrationType::stmfirstcolumn)
        ComputeShortTermVanillaPrice(isFreezeWeights,isInitMRS,isUpdateStrike);
	else if(itsMRSCalibType == ARM_MRSCalibrationType::diagstartfwd)
        ComputeDiagFWDVanillaPrice(isFreezeWeights,isInitMRS);

    /// Update caplet/floorlet prices
    ComputeImpliedCapletFloorletPrices(isFreezeWeights);

	if(itsIsFrontier)
		CreateAndSetCalibration_Frontier();

	if(itsModelType == ARM_PricingModelType::QGM1F && itsSkewCalFlag && itsSkewReCalFlag){
		/// To set the initial sigma Curve to avoid hedge jumb
		/// FIX FIX not correct if calibmethod has another modelparam at 0;
		/// please don't forget to change it as soon as possible
		GetPricingModel()->GetModelParams()->SetModelParam(GetOSWCalibMethod()->GetCalibParam(ARM_ModelParamType::Volatility));
	
        PreCalibrateSkew();
	}
}


////////////////////////////////////////////////////
///	Class   : ARM_CRFCalculator
///	Routines: Clone,toString
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_CRFCalculator::Clone() const
{
	return new ARM_CRFCalculator(*this);
}

string ARM_CRFCalculator::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream crfData;

    /// CRF Calculator specific datas viewing
    crfData <<indent <<"\n\n =======> CALLABLE REVERSE FLOATER CALCULATOR <====== \n\n";

    crfData << "\n"<<indent << "StartDate = " <<  itsStartDate.toString() << "\n";
    crfData << indent << "EndDate = " << itsEndDate.toString() << "\n";
    crfData << indent << "Pay/Rec = " << ARM_ParamView::GetMappingName(S_RECEIVE_PAY, itsPayRec) << "\n\n";

    crfData << indent << "Fixed Leg Datas :\n";
    crfData << indent << "EndDate = " <<  itsFixEndDate.toString() << "\n";
    crfData << indent << "Day Count	= " << ARM_ParamView::GetMappingName( S_DAYCOUNT,  itsFixDayCount ) << "\n\n";

    crfData << indent << "Reverse Leg Datas :\n";
    crfData << indent << "StartDate	= " <<  itsFixEndDate.toString() << "\n";
    crfData << indent << "Day Count	= " << ARM_ParamView::GetMappingName( S_DAYCOUNT,  itsCpnDayCount ) << "\n";
    crfData << indent << "Freq (Reset & Pay) = " << ARM_ParamView::GetMappingName(S_FREQUENCY, itsCpnFreq ) << "\n";
    crfData << indent << "Pay Calendar = " << itsCpnPayCal << "\n";
    crfData << indent << "Reset Timing = " << ARM_ParamView::GetMappingName(S_TIMING_MOD, itsCpnTiming ) << "\n";
	crfData << indent << "Stub Rule	= " << ARM_ParamView::GetMappingName(S_STUB_RULES, itsStubRule ) << "\n";
    crfData << indent << "Reset Gap	= " << itsCpnResetGap << "\n";
    crfData << indent << "Reset Calendar = " << itsCpnResetCal << "\n";
    crfData << indent << "Cpn Index Term = " << itsCpnIndexTerm << "\n";
    crfData << indent << "Cpn Index Day Count = " << ARM_ParamView::GetMappingName( S_DAYCOUNT, itsCpnIndexDayCount )  << "\n\n";

    crfData << indent << "Funding Leg Datas :\n";
    crfData << indent << "Day Count = " << ARM_ParamView::GetMappingName( S_DAYCOUNT,  itsFundDayCount ) << "\n";
    crfData << indent << "Frequency	= " << ARM_ParamView::GetMappingName(S_FREQUENCY, itsFundFreq ) << "\n\n";

    crfData << indent << "Exercise Datas :\n";
    crfData << indent << "Frequency (=Cpn freq) = " << ARM_ParamView::GetMappingName(S_FREQUENCY, itsCpnFreq ) << "\n";
    crfData << indent << "Notification Gap		= " << itsExerGap << "\n";
    crfData << indent << "Notification Calendar = " << itsExerCal << "\n";
    crfData << indent << "Non Call Period		= " << itsNbNonCall+(itsCpnTiming==K_ARREARS ? 1 : 0) << "\n\n";

    if(&GetGenSecurity())
    {
        const ARM_DealDescription dealDesc = GetGenSecurity()->GetDealDescription();
        size_t i,nbRows = dealDesc.GetRowsNb();

        crfData << indent << "Reverse Coupon Profile :\n";
        crfData << indent << "StartDate    \tStrike    \tLeverage  \tCpnMin    \tCpnMax    \tCapKAdj   \tFloorKAdj \n";
        for(i=1;i<nbRows;++i)
        {
            crfData << indent << ARM_Date(atof(dealDesc.GetElem(i,StartDate).c_str())).toString() << "\t";
            crfData << CC_NS(std,fixed) << atof(dealDesc.GetElem(i,Strike).c_str());
            if(dealDesc.GetElemFormat(i,Leverage) != ARM_MISSING_TYPE)
            {
                crfData << "\t" << CC_NS(std,fixed) << atof(dealDesc.GetElem(i,Leverage).c_str()) << "\t";
                crfData << CC_NS(std,fixed) << atof(dealDesc.GetElem(i,CpnMin).c_str()) << "\t";
                crfData << CC_NS(std,fixed) << atof(dealDesc.GetElem(i,CpnMax).c_str()) << "\t";
                crfData << CC_NS(std,fixed) << atof(dealDesc.GetElem(i,AdjCapSpread).c_str()) << "\t";
                crfData << CC_NS(std,fixed) << atof(dealDesc.GetElem(i,AdjFloorSpread).c_str()) << "\n";
            }
            else
                crfData << "\n";
        }

        crfData << "\n"<<indent << "Funding Spread Profile :\n";
        crfData << indent << "StartDate \tSpread   \n";
        for(i=1;i<nbRows;++i)
        {
            crfData << indent << ARM_Date(atof(dealDesc.GetElem(i,FundingStartDate).c_str())).toString() << "\t";
            crfData << CC_NS(std,fixed) << atof(dealDesc.GetElem(i,FundingSpread).c_str()) << "\n";
        }

        crfData << "\n"<< indent << "Nominal Profile :\n";
        crfData << indent << "PayDate   \tNominal   ";
	    if(itsFundNominal.GetDiscreteValues() != NULL && itsFundNominal.GetDiscreteValues()->size() > 0)
			crfData << "\tFunding Nominal   ";
		crfData << "\n";
		ARM_Date tmpDate;
        for(i=1;i<nbRows;++i)
        {
			tmpDate = ARM_Date(atof(dealDesc.GetElem(i,PayDate).c_str()));
            crfData << indent << tmpDate.toString() << "\t";
            crfData << CC_NS(std,fixed) << atof(dealDesc.GetElem(i,Nominal).c_str());
			if(itsFundNominal.GetDiscreteValues() != NULL && itsFundNominal.GetDiscreteValues()->size() > 0)
				crfData << "\t" << CC_NS(std,fixed) << const_cast< ARM_ReferenceValue& >(itsFundNominal).Interpolate(tmpDate.GetJulian());
			crfData << "\n";
        }

        crfData <<"\n"<< indent << "Notification Profile :\n";
        crfData << indent << "NotifDate \tFees      \n";
        for(i=1; i<nbRows; ++i)
        {
            crfData << indent << ARM_Date(atof(dealDesc.GetElem(i,EventDate).c_str())).toString() << "\t";
            crfData << CC_NS(std,fixed) << atof(dealDesc.GetElem(i,Fee).c_str()) << "\n";
		}
    }

	crfData << "\n\nCommon Calculator Part\n" << ARM_GenCalculator::toString(indent,nextIndent) << "\n\n";
    return crfData.str();
}


////////////////////////////////////////////////////
///	Class   : ARM_CRFCalculator
///	Routines: View
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////

void ARM_CRFCalculator::View(char* id, FILE* ficOut) const
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

    string indent = "";

    fprintf(fOut, "%s", toString(indent,"").c_str() );

        
    CC_Ostringstream crfData;
    crfData <<"\n\n"<< indent << "======> Auto Calibration Data <======\n\n";

    crfData <<"\n"<< indent << "Diagonal swaption portfolio to calibrate sigma Curve\n";
    crfData << indent << "Calibration=";
    crfData << indent << (itsOSWCalibFlag != ARM_SigmaCalibrationType::Unknown ? "On" : "Off");

    fprintf(fOut, "%s", crfData.str().c_str() );

    if(itsOSWCalibFlag != ARM_SigmaCalibrationType::Unknown)
    {
        ARM_StdPortfolioPtr diagonalSwaptionPF= GetOSWPortfolio();
        if(diagonalSwaptionPF!=ARM_StdPortfolioPtr(NULL))
            diagonalSwaptionPF->View(id,fOut);
    }
	CC_Ostringstream crfData0;
	crfData0 <<"\n"<< indent << "Diagonal swaption portfolio to calibrate Skew Curve\n\n";
	crfData << indent << "Calibration=";
     crfData << indent << (ARM_PricingModelType::QGM1F && itsSkewCalFlag ? "On" : "Off");

    fprintf(fOut, "%s", crfData0.str().c_str() );
    if(itsModelType == ARM_PricingModelType::QGM1F && itsSkewCalFlag)
    {
        if(itsOWSSkewPF)
            itsOWSSkewPF->View(id,fOut);
    }

    CC_Ostringstream crfData1;

    crfData1 <<"\n"<< indent << "Short term vanilla portfolio to calibrate MRS paramter\n";
	crfData1 << indent << "MRS Calib Type \t: " << ARM_ArgConvReverse_MRSCalibType.GetString( itsMRSCalibType ) << "\n";


    fprintf(fOut, "%s", crfData1.str().c_str() );

    if(itsMRSCalibType != ARM_MRSCalibrationType::Unknown)
    {
        ARM_StdPortfolioPtr vanillaPF= GetSTMPortfolio();
        if(vanillaPF!=ARM_StdPortfolioPtr(NULL))
            vanillaPF->View(id,fOut);
    }
    
    CC_Ostringstream crfData2;

    crfData2 <<"\n"<< indent << "Implied caplet & floorlet portfolio to calibrate strike spread adjument\n";
    crfData2 << indent << " Cap calibration=";
    crfData2 << indent << (itsCapCalibFlag ? "On" : "Off");
    crfData2 << indent << " Floor calibration=";
    crfData2 << indent << (itsFloorCalibFlag ? "On" : "Off");

    fprintf(fOut,"%s\n",crfData2.str().c_str());

    if(itsCapCalibFlag || itsFloorCalibFlag)
    {
        if(itsCapFloorPF)
            itsCapFloorPF->View(id,fOut);
    }

    CC_Ostringstream crfData3;
    crfData3 << indent << "Product to be priced = ";

    if(IsCRFToPrice())          crfData3 << indent << "Bermuda Reverse Floater";
    else if(IsCapToPrice())     crfData3 << indent << "Implicit Cap";
    else if(IsFloorToPrice())   crfData3 << indent << "Implicit Floor";
    else if(IsFundingToPrice()) crfData3 << indent << "Funding Leg";
    else if(IsStdLegToPrice())  crfData3 << indent << "Standard Coupon Leg (no min, no max cpn)";
    else if(IsRFLegToPrice())   crfData3 << indent << "Reverse Floater Coupon Leg";
    else if(IsStdSwapToPrice()) crfData3 << indent << "Underlying Standard Swap";
    else if(IsRFSwapToPrice())  crfData3 << indent << "Underlying Reverse Floater Swap";
    else if(IsBermudaToPrice()) crfData3 << indent << "Bermuda Standard Swaption";
    else crfData3 << indent << "Unknown !";
    fprintf(fOut,"\n\n %s\n",crfData3.str().c_str());

	//// ARRRRRGGGGGGGGG forced to const cast
    if(		itsConvexityModel != ARM_PricingModelPtr(NULL) 
		&&  const_cast<ARM_CRFCalculator*>(this)->itsConvexityModel != GetPricingModel() )
    {
        fprintf(fOut,"\n Convexity Model for in-arrears equivalent strike computation :\n");
        itsConvexityModel->View(id,fOut);
    }

    if ( ficOut == NULL )
       fclose(fOut);
}


double ARM_CRFCalculator::InterpolMeanRevParam(ARM_ReferenceValue* mrsParam)
{
	//// FIXME, we have to calculate MRS from the first start date
	//// after the first NotifDate
	double	res = MRS_DEFAULT_VALUE;

	if(mrsParam)
	{
		int oldExtrapolMeth = mrsParam->GetExtrapolMeth();
		mrsParam->SetExtrapolMeth(1); //forced to constant

		double	matu(0.0);
		ARM_Date asOf = GetMktDataManager()->GetAsOfDate();

		if (itsStartDate > asOf)
			matu = (itsEndDate-itsStartDate)/365.;
		else
			matu = (itsEndDate-asOf)/365.;

		res = mrsParam->CptReferenceValue(matu);
				
		mrsParam->SetExtrapolMeth(oldExtrapolMeth); // reset original method
	}

	return res;
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

