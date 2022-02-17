//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : CreditMetricsABSCDO.cpp
//
//   Description : Credit Metrics applied to ABS CDO
//
//   Date        : Jan 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CreditAsset.hpp"
#include "edginc/CreditMetricsABSCDO.hpp"
#include "edginc/CreditMetricsABSCDOLossCalculator.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/TrancheLossCalculator.hpp"
#include "edginc/FixedTrancheLossCalculator.hpp"
//#include "edginc/ConvolutionProduct.hpp"
#include "edginc/CreditTrancheLossConfig.hpp"
#include "edginc/CounterPartyCredit.hpp"

DRLIB_BEGIN_NAMESPACE

CreditMetricsABSCDO::~CreditMetricsABSCDO(){}

void CreditMetricsABSCDO::createLossCalculators(
        const DateTimeArray&                timeline,           /* (I) */
        CreditTrancheLossConfigConstSP      tranche,            /* (I) */
        CounterPartyCreditConstSP           cpty,               /* (I) */
        const DateTime&                     maturity,           /* (I) */
        Control*                            control,            /* (I) */
        Results*                            results,            /* (I) */
        bool                                recoverNotional,    /* (I) */
        YieldCurveConstSP                   discount,           /* (I) */
        IFixedTrancheLossCalculatorConstSP& lossCalculator,               /* (O) */
        IFixedTrancheLossCalculatorConstSP& recoveredNotionalCalculator,  /* (O) */
        IFixedTrancheLossCalculatorConstSP& conditionalLossCalc,          /* (O) */
        IFixedTrancheLossCalculatorConstSP& conditionalRecNtnlCalc) const //(O)
{
    try {
        if (!recoverNotional)
            throw ModelException("For ABS CDO, recover notional must be on");

        double lowStrike, highStrike;
        tranche->getTrancheStrikes(lowStrike, highStrike);
        double lossUnitToUse = calculateLossUnit(tranche);

        ITrancheLossCalculator* ctor;
        if (!useLossSkew && !useDecSkew && !useLossDecSkew) { // No base correlation
            DoubleArray lossBCOverride, decBCOverride, lossDecBCOverride;
            ctor = new CreditMetricsABSCDOLossCalculator(
                                                timeline, 
                                                tranche,
                                                lossUnitToUse,
                                                cpty,
                                                useSaddlePoint,
                                                numSaddlePoints,
                                                lossMarketFactors,
                                                decretionMarketFactors,
                                                useFastConvolution(tranche),
                                                lossBCOverride,
                                                decBCOverride,
                                                lossDecBCOverride,
                                                true);
        } 
        else { // Base correlation
            ctor = new ABSCDOBaseCorrelationLossCalculator(
                                                timeline, 
                                                tranche, 
                                                lossUnitToUse,
                                                cpty,
                                                maturity,
                                                useSaddlePoint,
                                                numSaddlePoints,
                                                lossMarketFactors,
                                                decretionMarketFactors,
                                                useFastConvolution(tranche),
												useLossSkew,
												useDecSkew,
												useLossDecSkew,
                                                authoriseNegativeEL);
        }
        lossCalculator.reset(
            ITrancheLossCalculator::createFixedTrancheLossCalculator(
                lowStrike, highStrike, ITrancheLossCalculatorConstSP(ctor)));
        
        ProxyCalculator* ptor = 
            new ProxyCalculator(dynamic_cast<ISupportProxy*>(ctor)); // XXX dynamic_cast

        recoveredNotionalCalculator.reset(
            ITrancheLossCalculator::createFixedTrancheLossCalculator(
                lowStrike, highStrike, ITrancheLossCalculatorConstSP(ptor)));
    } 
    catch (exception& e){
        throw ModelException(e, "CreditMetricsModel::createLossCalculators");
    }
}

/** Invoked by the containing model before the instrument data is fetched
    ie before CInstrument::GetMarket is invoked. Overrides
    CreditMetricsModel in order to get london floor curve */
void CreditMetricsABSCDO::preFetchMarketData(Model*            model,
                                             MarketDataConstSP market){
    //call the parent method first, so that any index map is populated
    CreditMetricsModel::preFetchMarketData(model,market);
    //now fetch the correlation
}

/** Called immediately after object constructed */
void CreditMetricsABSCDO::validatePop2Object() {
    CreditMetricsModel::validatePop2Object();
    if (useSaddlePoint && (numSaddlePoints > lossMarketFactors))
        throw ModelException("Can not use saddle point! numSaddlePoints > lossMarketFactors");
    if (lossMarketFactors <= 0 || lossMarketFactors > 100)
        throw ModelException("loss market factor numbers not in [1..100]");
    if (decretionMarketFactors <= 0 || decretionMarketFactors > 10)
        throw ModelException("decretion market factor numbers not in [1..10]");
}

/** Invoked when Class is 'loaded' */
void CreditMetricsABSCDO::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(CreditMetricsABSCDO, clazz);
    SUPERCLASS(CreditMetricsModel);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(useSaddlePoint,
                 "use saddle point instead of convolution");
    FIELD(numSaddlePoints,
                 "number of saddle points to explicitly calculate");
    FIELD(lossMarketFactors,
                 "number of market factor for loss");
    FIELD(decretionMarketFactors,
                 "number of market factor for decretion");
    
	FIELD(useLossSkew,
		"TRUE=use correlation skew for loss, "
		"FALSE=use beta for loss "
		"[default is FALSE]");
	FIELD_MAKE_OPTIONAL(useLossSkew);

	FIELD(useDecSkew,
		"TRUE=use correlation skew for decretion, "
		"FALSE=use beta for decretion "
		"[default is FALSE]");
	FIELD_MAKE_OPTIONAL(useDecSkew);

	FIELD(useLossDecSkew,
		"TRUE=use correlation skew for loss-decretion, "
		"FALSE=use beta for loss-decretion "
		"[default is FALSE]");
	FIELD_MAKE_OPTIONAL(useLossDecSkew);

    FIELD(authoriseNegativeEL,
        "TRUE=authorise negative expected losses, "
        "FALSE=do not authorise negative expected losses (floors to zero) "
        "[default is FALSE]");
    FIELD_MAKE_OPTIONAL(authoriseNegativeEL);
}

/** Private constructor */
CreditMetricsABSCDO::CreditMetricsABSCDO(const CClassConstSP& clazz) : 
    CreditMetricsModel(clazz), useSaddlePoint(true), 
        lossMarketFactors(0), decretionMarketFactors(0),
        useLossSkew(false), useDecSkew(false), useLossDecSkew(false), 
        authoriseNegativeEL(false) {}

/** Default constructor */
IObject* CreditMetricsABSCDO::defaultConstructor(){
    return new CreditMetricsABSCDO(CreditMetricsABSCDO::TYPE);
}

CClassConstSP const CreditMetricsABSCDO::TYPE = 
CClass::registerClassLoadMethod(
    "CreditMetricsABSCDO", typeid(CreditMetricsABSCDO), load);

bool CreditMetricsABSCDOLoad() {
    return (CreditMetricsABSCDO::TYPE != 0);
}
DRLIB_END_NAMESPACE
