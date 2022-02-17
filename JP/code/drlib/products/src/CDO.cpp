//----------------------------------------------------------------------------
//
// Group       : Derivatives Research
//
// Filename    : CDO.cpp
//
// Description : Collateralized Debt Obligation Instrument. An instance of
//               a GeneralisedCDO where the underlying portfolio is a tranche
//
// Author      : Mark A Robson/Sebastien Hitier/Antoine Gregoire
//
// Date        : 28 Oct 2004
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/CDO.hpp"
#include "edginc/CCMPriceUtil.hpp"
#include "edginc/SingleCreditAsset.hpp"
#include "edginc/OutputRequestUtil.hpp"
#include "edginc/CreditCashFlow.hpp"
#include "edginc/CounterPartyCredit.hpp"
#include "edginc/CtgLegLossPerDefault.hpp"
#include "edginc/FeeLegReductionPerDefault.hpp"
#include "edginc/CreditTrancheLossConfig.hpp"
#include "edginc/CreditContingentLegBase.hpp"


DRLIB_BEGIN_NAMESPACE

CDO::~CDO()
{}

CDO::CDO(CClassConstSP clazz) :
	GeneralisedCDO(clazz)
{}

/** Public constructor */
CDO::CDO(double                       lowStrike,
         double                       highStrike,
         bool                         isLong,
         CBoolSP                      recoverNotional,
         smartPtr<CounterPartyCredit> cptyInfo,
         ICreditContingentLegSP       cLeg,
         ICreditFeeLegSP              fLeg,
         CDOPortfolioSP               portfolio,
         const string&                yieldCurveName,
         CIntSP                       triggerDelay,
         CIntSP                       defaultToCalculationDelay,
         DateTime                     lastTriggerDate,
         CDoubleSP                    temporaryLossAmount,
         const string&                settlementBDC,
         HolidayWrapper               settlementHols) :
    GeneralisedCDO(TYPE,
                   isLong,
                   cptyInfo,
                   cLeg,
                   fLeg,
                   portfolio,
                   recoverNotional,
                   yieldCurveName,
                   triggerDelay,
                   defaultToCalculationDelay,
                   lastTriggerDate,
                   temporaryLossAmount,
                   settlementBDC,
                   settlementHols),
    lowStrike(lowStrike),
    highStrike(highStrike)
{
    validatePop2Object();
}

void CDO::validatePop2Object() {
    static const string method("CDO::validatePop2Object");

    CreditTrancheLossConfig* portf =
        dynamic_cast<CreditTrancheLossConfig*>(portfolio.get());

    if (!portf) {
        // Got something other than a CreditTrancheLossConfig - replace it
        // with a CreditTrancheLossConfig
        portfolio = ICreditLossConfigSP(new CreditTrancheLossConfig(
            portfolio->getName(),
            getValueDate(),
            lowStrike,
            highStrike,
            portfolio));
    }

    GeneralisedCDO::validatePop2Object();
}


IObject* CDO::defaultConstructor() {
    return new CDO();
}


/** Invoked once at start up when this class is 'loaded' */
void CDO::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(CDO, clazz);
    SUPERCLASS(GeneralisedCDO);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(lowStrike,  "Low strike (absolute value, NOT a percentage)");
    FIELD(highStrike, "High strike (absolute value, NOT a percentage)");
}

CClassConstSP const CDO::TYPE = CClass::registerClassLoadMethod(
    "CDO", typeid(CDO), load);



//------------------------------------------------------------------------------
//                                CDO::QuickPricer
//------------------------------------------------------------------------------
/** Invoked once at start up when this class is 'loaded' */
void CDO::QuickPricer::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(QuickPricer, clazz);
    SUPERCLASS(KComponent);
    IMPLEMENTS(ConvolutionEngine::IIntoProduct);
	IMPLEMENTS(FDModel::IIntoProduct);
    IMPLEMENTS(Theta::Shift);
    IMPLEMENTS(LastSensDate);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD           (today, "value date");
    FIELD_MAKE_OPTIONAL    (today);
    FIELD           (tradeDate, "Trade start date");
    FIELD_MAKE_OPTIONAL(tradeDate);
    FIELD_MAKE_NONTWEAKABLE(today);
    FIELD           (lowStrike, "Low strike (absolute value, NOT a percentage)");
    FIELD_MAKE_NONTWEAKABLE(lowStrike);
    FIELD           (highStrike, "High strike (absolute value, NOT a percentage)");
    FIELD_MAKE_NONTWEAKABLE(highStrike);
    FIELD           (isLong, "Long or Short");
    FIELD_MAKE_NONTWEAKABLE(isLong);
    FIELD                  (recoverNotional, "Recover notional on the fee leg. Default: false");
    FIELD_MAKE_OPTIONAL    (recoverNotional);
    FIELD                  (cptyInfo, "cpty charge info");
    FIELD_MAKE_OPTIONAL    (cptyInfo);
    FIELD_MAKE_NONTWEAKABLE(cptyInfo);
    FIELD                  (portfolio, "tranche portfolio");
    FIELD_MAKE_NONTWEAKABLE(portfolio);
    FIELD                  (legConventions, "Conventions that define fee and contingent legs");
    FIELD_MAKE_NONTWEAKABLE(legConventions);
    FIELD                  (expiry, "Tranche expiry");
    FIELD_MAKE_NONTWEAKABLE(expiry);
    FIELD           (spread, "Spread (fee leg fixed coupon)");
    FIELD_MAKE_NONTWEAKABLE(spread);
    FIELD           (upfrontPayment, "Upfront fee payment (default is 0)");
    FIELD_MAKE_OPTIONAL    (upfrontPayment);
    FIELD           (legNotional, "Fee and contingent leg notional "
                            "(default is 1). The leg notional can be overriden "
                            "via feeLegNotionalOverride.");
    FIELD_MAKE_OPTIONAL    (legNotional);
    FIELD                  (feeLegNotionalOverride, "Fee leg notional override");
    FIELD_MAKE_OPTIONAL    (feeLegNotionalOverride);

    FIELD                  (triggerDelay, "Delay between market default and "
                            "eventDeterminationDate");
    FIELD_MAKE_OPTIONAL    (triggerDelay);
    FIELD                  (defaultToCalculationDelay, "Delay between market "
                            "default and the associated calculation date");
    FIELD_MAKE_OPTIONAL    (defaultToCalculationDelay);
    FIELD           (lastTriggerDate, "Last date when a default occurred "
                            "during the protection period can be triggered. "
                            "Default: Protection end date.");
    FIELD_MAKE_OPTIONAL    (lastTriggerDate);
    FIELD                  (temporaryLossAmount, "Assumption for the loss amount "
                            "associated to a default, used to determine the fees "
                            "in accrual periods between the credit event and the "
                            "calculation date. Default: 100%");
    FIELD_MAKE_OPTIONAL    (temporaryLossAmount);
    FIELD           (settlementBDC, "Bad day convention for the tranche. "
                            "Used to adjust the delays regarding defaults' "
                            "settlements. Default: None");
    FIELD_MAKE_OPTIONAL    (settlementBDC); // to be made mandatory ideally
    FIELD           (settlementHols, "Holidays for the tranche. Used to "
                            "adjust the delays regarding defaults' settlements. "
                            "Default: weekends only");
    FIELD_MAKE_OPTIONAL    (settlementHols); // to be made mandatory ideally

    FIELD_NO_DESC                     (cdo);
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(cdo);
}


IObject* CDO::QuickPricer::defaultConstructor() {
    return new QuickPricer();
}

/** Private constructor (for reflection) */
CDO::QuickPricer::QuickPricer() :
    KComponent("","GeneralisedCDO",TYPE),
    upfrontPayment(0.0),
    legNotional(1.0),
    settlementBDC("N")
{}

/** Called immediately after object constructed */
void CDO::QuickPricer::validatePop2Object() {
    try {
		discount = YieldCurveWrapper(legConventions->getDiscountName());
		
        // init cdo
        cdo.reset(new CDO(lowStrike,
                          highStrike,
                          isLong,
                          recoverNotional,
                          cptyInfo,
                          ICreditContingentLegSP(   ),
                          ICreditFeeLegSP(   ),
                          portfolio,
                          legConventions->getDiscountName(),
                          triggerDelay,
                          defaultToCalculationDelay,
                          lastTriggerDate,
                          temporaryLossAmount,
                          settlementBDC,
                          settlementHols));
    }
    catch (exception& e) {
        throw ModelException(e, "CDO::QuickPricer::validatePop2Object");
    }
}

/** Retrieve market data from cache for all our components */
void CDO::QuickPricer::GetMarket(const IModel*       model,
                                 const CMarketDataSP market) {
    try {

		// check if we are using KComponent framework
		//if(FDModel::TYPE->isInstance(model))
		if(false)
		{
			
			// Call getMarket on underlying KComponent - needed for tree
			// for the moment only call this when the pricing model is of tree type
			// could remove this later if feel more confident with KComponent
			KComponent::GetMarket(model, market);

		}
		else
		{
			market->GetReferenceDate(today);
			legConventions->getMarket(model, market.get());

            if(tradeDate.empty()) tradeDate = today;

			// build fee / contingent legs
			DateTime startDate = legConventions->startDate(tradeDate);
			DateTime maturity = expiry->toDate(tradeDate);

			// build fee leg
			cdo->fLeg = legConventions->generateFeeLeg(
					startDate, maturity, spread, upfrontPayment,
					(!feeLegNotionalOverride ?
					legNotional : // if override not present, use legNotional
					feeLegNotionalOverride->doubleValue()));

			// build contingent leg
			cdo->cLeg = legConventions->generateContingentLeg(
					startDate, maturity, spread, upfrontPayment, legNotional);

			cdo->GetMarket(model, market);
		}
    } catch (exception& e) {
        throw ModelException(e, "CDO::QuickPricer::GetMarket");
    }
}

// Setup the component. 
// The derived function must also call KComponent::setup() and the setup() of its underlyings
void CDO::QuickPricer::setup(const IModel* model, const MarketData* market)
{
	static const string method = "CDO::QuickPricer::setup";
	
	KComponent::setup(model, market);
	
	market->GetReferenceDate(today);
	legConventions->getMarket(model, market);
	
     if(tradeDate.empty()) tradeDate = today;

	// build fee / contingent legs
	DateTime startDate = legConventions->startDate(tradeDate);
	DateTime maturity = expiry->toDate(tradeDate);

   

	// build fee leg
	cdo->fLeg = legConventions->generateFeeLeg(
		startDate, maturity, spread, upfrontPayment,
		(!feeLegNotionalOverride ?
		legNotional : // if override not present, use legNotional
		feeLegNotionalOverride->doubleValue()));

	// build contingent leg
	cdo->cLeg = legConventions->generateContingentLeg(
		startDate, maturity, spread, upfrontPayment, legNotional);

	// cascade the setup call
	cdo->setup(model, market);

}

/** Called after market data has been retrieved */
void CDO::QuickPricer::Validate() {
    cdo->Validate();
}

//// Required part of CInstrument
DateTime CDO::QuickPricer::getValueDate() const {
    return cdo->getValueDate();
}


//// CDO is guaranted to have a CDOPortfolio inside
CDOPortfolioSP CDO::getCDOPortfolio() const
{
    // cast the CreditLossConfig into a CreditTranchelossConfig
    CreditTrancheLossConfig *tlc    = DYNAMIC_CAST(CreditTrancheLossConfig, portfolio.get());
    ICreditLossConfigConstSP ptf      = tlc->getInnerLossConfig();
    return CDOPortfolioSP( DYNAMIC_CAST(CDOPortfolio,ptf.get()));
}

//// function used only in the TrancheIndexLeastSquareFit
void CDO::getMaturityAndStrikesPercent(DateTime& t, double &k1, double &k2) const
{
    // cast the CreditLossConfig into a CreditTranchelossConfig
    CreditTrancheLossConfig *tlc = DYNAMIC_CAST(CreditTrancheLossConfig, portfolio.get());

    // get the portfolio long notional from the CreditTrancheLossConfig
    double portfolioLongNotional = tlc->portfolioNotional();

    k1 = lowStrike  / portfolioLongNotional;
    k2 = highStrike / portfolioLongNotional;
    t  = lastObservationDate();
}


/** Returns the name of the instrument's discount currency. */
string CDO::QuickPricer::discountYieldCurveName() const {
    return cdo->discountYieldCurveName();
}

/**
 * when to stop tweaking (need to change infrastructure to route through
 * model rather than instrument)
 * */
DateTime CDO::QuickPricer::endDate(const Sensitivity* sensControl) const {
    return cdo->endDate(sensControl);
}

//// Required part of Theta::Shift
bool CDO::QuickPricer::sensShift(Theta* shift) {
    return cdo->sensShift(shift);
}

/** Creates an instance of an ConvolutionProduct */
IGeneralisedConvolutionProduct* CDO::QuickPricer::createProduct(
    ConvolutionEngineConstSP model) const
{
    return cdo->createProduct(model);
}

/** Creates an instance of an FDProduct (tree product) */
FDProductSP CDO::QuickPricer::createProduct( FDModel * model ) const
{
	return cdo->createProduct(model);
}

CClassConstSP const CDO::QuickPricer::TYPE =
    CClass::registerClassLoadMethod("CDO::QuickPricer",
                                    typeid(CDO::QuickPricer),
                                    load);

/* For class loading (avoid having header file) */
bool CDOLoad() {
    return (CDO::TYPE != 0);
}

DRLIB_END_NAMESPACE
