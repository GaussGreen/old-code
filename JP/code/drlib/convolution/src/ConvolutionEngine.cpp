//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Description : Family of CModels that use convolution to get an 
//                 'effective curve'
//
//   Date        : 18th Nov 2005
//
//   Author      : Mark Robson
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SingleCreditAsset.hpp"
#include "edginc/ConvolutionEngine.hpp"
#include "edginc/ConvolutionModel.hpp"
#include "edginc/ConvolutionProduct.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/Control.hpp"
#include "edginc/MarketDataFetcher.hpp"
#include "edginc/SRMFXDiffuse.hpp"
#include "edginc/QuantoCDSAlgorithm.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/Addin.hpp"
#include "edginc/RhoParallel.hpp"
#include "edginc/RhoPointwise.hpp"
#include "edginc/ICreditLossConfig.hpp"
#include "edginc/IEffectiveCurveLossGen.hpp"
#include "edginc/CreditTrancheLossConfig.hpp"
#include "edginc/FlatCDO2LossConfig.hpp"
#include "edginc/IModelConfigMapper.hpp"
#include "edginc/TrancheLossCurveGen.hpp"
#include "edginc/IForwardRatePricer.hpp"
#include "edginc/FtDEffectiveCurveGen.hpp"

DRLIB_BEGIN_NAMESPACE

/** Key method providing access to the pricer */
IForwardRatePricerSP ConvolutionEngine::getForwardRatePricer() const
{
    static const string method = "ConvolutionEngine::getForwardRatePricer";

    // retrieve the model to be used in the calculation
    // of fee leg forward rates (from this model)
    IHasForwardRatePricer* ihfrp = dynamic_cast<IHasForwardRatePricer*>
        (const_cast<IConvolutionModel*>(convolutionModel.get()));
    if (!ihfrp)
    {
        throw ModelException(method,
            "Model must implement IHasForwardRateModel");
    }
    IForwardRatePricerSP frModel = ihfrp->getForwardRatePricer();

    return frModel;
}

ConvolutionEngine::IIntoProduct::~IIntoProduct()
{}

ConvolutionEngine::~ConvolutionEngine()
{}

/** overridden to [reference] copy QuantoCDSAlgorithmSP/irGridPtsCache */
IObject* ConvolutionEngine::clone() const{
    ConvolutionEngine* myCopy = DYNAMIC_CAST(ConvolutionEngine, Model::clone());
    myCopy->quantoAlgorithm = quantoAlgorithm; // just copy reference
    myCopy->irGridPtsCache = irGridPtsCache; // just copy reference
    return myCopy;
}

//// called immediately after object constructed
void ConvolutionEngine::validatePop2Object(){
    static const string method = "ConvolutionEngine::validatePop2Object";
    if (creditChargeViewType != ConvolutionProduct::CC_CONSERVATIVE &&
        creditChargeViewType != ConvolutionProduct::CC_AGGRESSIVE){
        throw ModelException(method, "Counterparty: creditCharge View "
                             "type flag must be "+
                             ConvolutionProduct::CC_CONSERVATIVE+" or "+ 
                             ConvolutionProduct::CC_AGGRESSIVE);
    }
}

/** Overrides default Model implementation to retrieve today, then
    calls ConvolutionModel::preFetchMarketData before invoking
    parent method, and then finally invoking 
    ConvolutionModel::postFetchMarketData */
void ConvolutionEngine::getInstrumentsAndModelMarket(
    MarketDataConstSP       market,
    IInstrumentCollectionSP insts)
{
    static const string method("ConvolutionEngine::"
                               "getInstrumentsAndModelMarket");
    try {
        fetchToday(market); // populate today in model
        // allow ConvolutionModel early access to market data
        convolutionModel->preFetchMarketData(this, market); 
        // call method that we're overriding
        Model::getInstrumentsAndModelMarket(market, insts);
        // and in case ConvolutionModel has identifed extra data it needs
        convolutionModel->postFetchMarketData(this, market); 
    } catch (exception& e){
        throw ModelException(e, "ConvolutionEngine::"
                             "getInstrumentsAndModelMarket");
    }
}


/** calculate single price and store result in CResult */
void ConvolutionEngine::Price(CInstrument*         instrument, 
                              CControl*            control, 
                              CResults*            results)
{
    static const string method("ConvolutionEngine::Price");

    // retrieve the model to be used in the calculation
    // of fee leg forward rates (from this model)
    IForwardRatePricerSP frModel = getForwardRatePricer();

    auto_ptr<IGeneralisedConvolutionProduct> prod(createProduct(instrument));
    bool isPricing = control->isPricing();
    if (isPricing && quantoAlgorithm.get()) {
        // tell quanto algorithm to save where it's interpolating the vol
        // (Used for IR vega tweaks)
        irGridPtsCache.reset(new IRGridPointCache(true));
        quantoAlgorithm->setIRGridPointCache(irGridPtsCache);
    }

    prod->price(convolutionModel.get(), control, results, frModel);

    if (isPricing && quantoAlgorithm.get()) {
        // switch off caching of where we've interpolated
        irGridPtsCache->setCachingMode(false);
    }
}

/** Creates IGeneralisedConvolutionProduct and returns it */
IGeneralisedConvolutionProduct* ConvolutionEngine::createProduct(
    const CInstrument* instrument) const
{
    static const string method("ConvolutionEngine::createProduct");
    try {
        if (!IIntoProduct::TYPE->isInstance(instrument)){
            throw ModelException(method, "Instrument of type "+
                                 instrument->getClass()->getName() + " "
                                 "does not support being priced with "
                                 "ConvolutionEngine");
        }
        const IIntoProduct& intoProd = 
            dynamic_cast<const IIntoProduct&>(*instrument);
        return intoProd.createProduct(ConvolutionEngineConstSP(this));
    } 
    catch (exception& e) {
        throw ModelException(e, method, "Failed to create product");
    }
}


/** Creates a IEffectiveCurveLossGen for the supplied ICreditLossConfig object. 
    The ICreditLossConfig object must implement the IntoLossGen interface */
// NB: "mapper" is currently NOT USED
ICreditLossGenSP ConvolutionEngine::lossGenerator(
    ICreditLossConfigConstSP lossConfig,
    IModelConfigMapperConstSP mapper) const
{
    static const string method("ConvolutionEngine::lossGenerator");

    // very simple implementation here although more sophisticated ones
    // could pass a different EffectiveLossCurveModelConfig depending
    // on the lossConfig's name.
    const IIntoLossGen* intoLossGen = 
        dynamic_cast<const IIntoLossGen*>(lossConfig.get());

    if (!intoLossGen){
        throw ModelException(method, "The instrument's loss config '" +
                             lossConfig->getName() +
                             "' can not be priced using this model.");
    }
    
    return intoLossGen->lossGenerator(ConvolutionEngineConstSP::attachToRef(this));
}


/** Creates an IEffectiveCurveGen for the supplied ICreditLossConfig object. 
    The ICreditLossConfig object must implement the IntoEffCurveGen interface */
// NB: "mapper" is currently NOT USED
ICreditLossGenSP ConvolutionEngine::effCurveGenerator(
    ICreditLossConfigConstSP  lossConfig,
    CounterPartyCreditConstSP cpty,
    const bool                recoverNotional,
    IModelConfigMapperConstSP mapper) const
{
    static const string method("ConvolutionEngine::effCurveGenerator");

    // very simple implementation here although more sophisticated ones
    // could pass a different EffectiveLossCurveModelConfig depending
    // on the lossConfig's name.
    const IIntoEffCurveGen* intoEffCurveGen = 
        dynamic_cast<const IIntoEffCurveGen*>(lossConfig.get());

    if (!intoEffCurveGen){
        throw ModelException(method, "The instrument's loss config '" +
                             lossConfig->getName() +
                             "' can not be priced using this model.");
    }
    
    return intoEffCurveGen->effCurveGenerator(
        ConvolutionEngineConstSP::attachToRef(this),
        cpty,
        recoverNotional,
        mapper);
}


/** Creates an IEffectiveCurveLossGen for a 'tranche' */
IEffectiveCurveLossGenSP ConvolutionEngine::createLossGenerator(
    CreditTrancheLossConfigConstSP trancheLossCfg) const
{
    return TrancheLossCurveGenSP(new TrancheLossCurveGen(
        convolutionModel, 
        ConvolutionEngineConstSP::attachToRef(this), 
        trancheLossCfg));
}

/** Creates an IEffectiveCurveLossGen for a 'tranche' */
IEffectiveCurveLossGenSP ConvolutionEngine::createLossGenerator(
    FlatCDO2LossConfigConstSP tsLossCfg) const
{
    return FlatCDO2LossCurveGenSP(new FlatCDO2LossCurveGen(
        convolutionModel, 
        ConvolutionEngineConstSP::attachToRef(this), 
        tsLossCfg));
}

/** Creates an IEffectiveCurveGen for an 'NtD' */
ICreditLossGenSP ConvolutionEngine::createEffCurveGenerator(
    NToDefaultLossConfigConstSP ntdLossCfg,
    CounterPartyCreditConstSP   cpty,
    const bool                  recoverNotional) const
{
    return convolutionModel->createEffCurveGenerator(ntdLossCfg,
                                                     cpty,
                                                     recoverNotional);
}


/** Whether risk mapping is applicable/enabled */
IModel::WantsRiskMapping ConvolutionEngine::wantsRiskMapping() const {
    return convolutionModel->wantsRiskMapping();
}

ConvolutionEngine::ConvolutionEngine(): 
    CModel(TYPE), doPVToSpot(false),
    // sensible default for ir vol interpolation (only used for
    // quanto adjustment)
    calibrationStyle("CMS"), calibrationMaturity("10Y")
{}

/** Constructor - mainly for backward compatibility. Allow models to create
    this object and then run the pricing through this class */
ConvolutionEngine::ConvolutionEngine(IConvolutionModelSP  convolutionModel,
                                     const DateTime&      today,
                                     const string&        creditChargeViewType,
                                     bool                 pvToSpot,  
                                     const DateTime&      cfCutOffDate,
                                     const string& calibrationStyle,
                                     const string& calibrationMaturity):
    CModel(TYPE, today), convolutionModel(convolutionModel),
    creditChargeViewType(creditChargeViewType),
    doPVToSpot(pvToSpot), theCfCutOffDate(cfCutOffDate),
    calibrationStyle(calibrationStyle),
    calibrationMaturity(calibrationMaturity),
    quantoAlgorithm(quantoAlgorithm) {
    validatePop2Object();
}

/** Overridden and redirected to ConvolutionModel */
MarketDataFetcherSP ConvolutionEngine::createMDF() const{
    return convolutionModel->createMDF();
}

/** Returns all the points on the skew surface to which this
    model for the supplied instrument is sensitive.
    A null return value is ok - it is interpreted as the greek being 
    NotApplicable */
BetaSkewGridPointArrayConstSP ConvolutionEngine::getSensitiveBetaSkewPoints(
    OutputNameConstSP  outputName,
    const CInstrument* inst,
    bool               tweakAll, 
    const int          numberOfNeighbours) const
{
    static const string method("ConvolutionEngine::getSensitiveBetaSkewPoints");

    auto_ptr<IGeneralisedConvolutionProduct> product(createProduct(inst));

    ICreditLossConfigConstSP lossConfig = product->getLossConfig();
    ICreditLossGenSP lossGen = lossGenerator(lossConfig,
                                             IModelConfigMapperConstSP());
    IEffectiveCurveLossGenSP effCurveLossGen =
        DYNAMIC_POINTER_CAST<IEffectiveCurveLossGen>(lossGen);

    if (!effCurveLossGen) {
        throw ModelException(method,
                             "Internal error: The loss generator is not of "
                             "type IEffectiveCurveLossGen");
    }
    const DateTime& maturity = product->lastObservationDate();
    return effCurveLossGen->getSensitiveBetaSkewPoints(outputName, 
                                                       maturity, 
                                                       product->getDiscount(),
                                                       tweakAll,
                                                       numberOfNeighbours);
}

/** Returns conservative/aggressive etc. To do: either enums or sort out
    where strings defined */
const string& ConvolutionEngine::getCreditChargeViewType() const {
    return creditChargeViewType;
}

/** Kapital backward compatibility mode - to be removed eventually.
    Returns whether we pv to cfCutOffDate() */
bool ConvolutionEngine::pvToSpot() const{
    return doPVToSpot;
}

/** Kapital backward compatibility mode - to be removed eventually.
    Ignore cashflows as well as protection etc before this date. */
const DateTime& ConvolutionEngine::cfCutOffDate() const{
    return theCfCutOffDate;
}

/** Returns true if the Sensitivity shifts yield curve type parameters */
bool ConvolutionEngine::isYCSens(const Sensitivity* sensitivity){
    // implementation is poor - should make IShifts derive from common class
    return (RhoParallel::TYPE->isInstance(sensitivity) ||
            RhoPointwise::TYPE->isInstance(sensitivity));
}

/** Returns true if the Sensitivity shifts IR Vol Params */
bool ConvolutionEngine::isIRVegaPointwiseSens(const Sensitivity* sensitivity) {
    // implementation is poor - should make IShifts derive from common class
    return (IRVegaPointwise::TYPE->isInstance(sensitivity));
}

/** Returns the end date for rho tweaks (as far as the quanto algorithm
    is concerned) - needs tidying up */
DateTime ConvolutionEngine::rhoEndDate(const Sensitivity* sensitivity,
                                       const DateTime&    irVegaEndDate) const
{
    // we need a better way of ensuring that we cover any other similar
    // sensitivities to rho pointwise
    if (quantoAlgorithm.get() && RhoPointwise::TYPE->isInstance(sensitivity)) {
        // only for quanto case
        OutputNameConstSP outputName =
            dynamic_cast<const RhoPointwise*>(sensitivity)->getMarketDataName();
        return irGridPtsCache->rhoEndDate(outputName, irVegaEndDate);
    }
    return irVegaEndDate;
}

/** when to stop tweaking for Yield Curve type tweaks */
DateTime ConvolutionEngine::endDateForYCTweaks(
    const IGeneralisedConvolutionProduct* product,
    IEffectiveCurveLossGenSP              effCurveLossGen,
    const Sensitivity*                    sensitivity,
    const DateTime&                       creditEndDate) const
{
    DateTime maxDate(getValueDate());
    maxDate = product->lastYCSensDate(maxDate);

    // check expiries for each name
    maxDate = maxDate.max(effCurveLossGen->maxLossConfigsTweakableDate(creditEndDate));

    // see if quanto goes further out
    return rhoEndDate(sensitivity, maxDate);
}

/** When to stop tweaking - to do change infrastructure to route through
    model rather than instrument */
DateTime ConvolutionEngine::endDate(const Instrument* instrument, 
                                    const Sensitivity* sensControl) const 
{
    static const string method("ConvolutionEngine::endDate");

    // to do: needs sorting out. In particular, redo when new quanto algorithm
    // is implemented.

    // create product
    auto_ptr<IGeneralisedConvolutionProduct> product(createProduct(instrument));

    ICreditLossConfigConstSP lossConfig = product->getLossConfig();
    ICreditLossGenSP lossGen = lossGenerator(lossConfig,
                                             IModelConfigMapperConstSP());
    IEffectiveCurveLossGenSP effCurveLossGen =
        DYNAMIC_POINTER_CAST<IEffectiveCurveLossGen>(lossGen);
    if (!effCurveLossGen) {
        throw ModelException(method,
                             "Internal error: The loss generator is not of "
                             "type IEffectiveCurveLossGen");
    }
    const DateTime& lastObservationDate = product->lastObservationDate();
    DateTimeArraySP timeline = effCurveLossGen->generateTimeline(lastObservationDate);
    const DateTime creditEndDate = timeline->back();

    // check sensControl type - needs tidying up
    if (dynamic_cast<const CreditTweak *>(sensControl)) {
        return creditEndDate;
    }

    const DateTime& dateForYCTweaks = endDateForYCTweaks(
        product.get(), effCurveLossGen, sensControl, creditEndDate);

    if (isYCSens(sensControl)) {
        return dateForYCTweaks;
    }
   
    // be safe
    DateTime maxDate(dateForYCTweaks);
    maxDate = creditEndDate.max(maxDate); // Redundant?
    if (isIRVegaPointwiseSens(sensControl)){
        // there's some weird with the quanto algorithm - the dependency
        // seems to go beyond what you'd expect. Very hack solution - move
        // forward one year
        return MaturityPeriod::toDate(1, "A", maxDate);
    }
    return maxDate;
}


/** Returns an instance of IAlgorithm. Typically the quantoCDSParSpreads
    parameter would be ignored but is there in case you want to switch
    the algorithm dependent upon some property of the quanto'd curve */
QuantoCDSParSpreads::IAlgorithmSP ConvolutionEngine::cdsQuantoAlgorithm(
    const QuantoCDSParSpreads* quantoCDSParSpreads) const{
    if (!quantoAlgorithm){
        createQuantoAlgorithm(false);
    }
    return quantoAlgorithm;
}
    
/** sets debug state to specified value. If true, then any calls to
    cdsQuantoAlgorithm will return an object that will cache debug data.
    This can be retrieved via getDebugInfo() */
void ConvolutionEngine::selectDebugState(bool switchOn){
    createQuantoAlgorithm(switchOn);
}

/** Returns debug info - may be null eg if the algorithm has not been 
    used */
IObjectSP ConvolutionEngine::getDebugInfo() const{
    if (!quantoAlgorithm){
        return IObjectSP();
    }
    return quantoAlgorithm->getDebugInfo();
}

/** Builds QuantoCDSAlgorithm object */
void ConvolutionEngine::createQuantoAlgorithm(bool debugOn) const {
    // decent error message in case empty strings have been passed in
    if (calibrationStyle.empty() || calibrationMaturity.empty()){
        throw ModelException("ConvolutionEngine:createQuantoAlgorithm", "Both "
                             "calibrationStyle and calibrationMaturity model"
                             " parameters must be specified when cds curves "
                             "need to be quanto'd");
    }
    /* Note: the last 4 strings below are for the 'correlation swap'. The
       actual values seem to matter very little. The day count convention and
       payment frequence match maxim/aladdin. For the expiry/maturity
       maxim/aladdin uses 10Y into 10Y. Switching to that caused no numbers to
       change except that several tests failed which had very very large spot
       vols - so I left it as it is here */
    quantoAlgorithm.reset(new QuantoCDSAlgorithm("1F-STANDARD", // hard coded
                                                 calibrationStyle,
                                                 calibrationMaturity,
                                                 false, // skipIRBadVols
                                                 SRMFXDiffuse::USE_LAST_LEVEL,
                                                 0.0, // fx cut off level
                                                 false, // floor negative vols
                                                 // swap starting in 1Y
                                                 // maturity 5Y
                                                 "1Y", "5Y", "Act/365F", "A",
                                                 debugOn));
}

/** Uses quanto calibration parameters to identify relevant points */
IRGridPointAbsArraySP ConvolutionEngine::getSensitiveIRVolPoints(
    OutputNameConstSP  outputName,
    const CInstrument* inst) const{
    if (!irGridPtsCache){
        // this can happen if eg an IRVol is directly embedded into a Yield
        // Curve (rather than being in the market data cache)
        // Return empty array rather than fail
        return IRGridPointAbsArraySP(new IRGridPointAbsArray());
    }
    return irGridPtsCache->sensitiveIRVolPoints(outputName);
}


/** overrides Models implementation to shift cfCutOffDate */
bool ConvolutionEngine::sensShift(Theta* shift){
    // call base method first
    CModel::sensShift(shift);
    if (!theCfCutOffDate.empty()){
        // alter immediate data
        theCfCutOffDate = shift->rollDate(theCfCutOffDate);
    }
    return true; // then shift components
}

IObject* ConvolutionEngine::defaultConstructor(){
    return new ConvolutionEngine();
}

void ConvolutionEngine::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(ConvolutionEngine, clazz);
    SUPERCLASS(CModel);
    IMPLEMENTS(QuantoCDSParSpreads::IAlgorithmBuilder);
    IMPLEMENTS(IEffectiveCurveLossModelConfig);
    // JLHP why are the other interfaces implemented here not registered?
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(convolutionModel, "drives way probabilities are computed "
          "in convolution");
    FIELD(creditChargeViewType, "credit charge view type, either "
                 "CONSERVATIVE or AGGRESSIVE");
    FIELD_NO_DESC(doPVToSpot);
    FIELD_MAKE_TRANSIENT(doPVToSpot);
    FIELD_NO_DESC(theCfCutOffDate);
    FIELD_MAKE_TRANSIENT(theCfCutOffDate);
    FIELD(calibrationStyle, "calibration style (eg CMS) for quanto "
                 "adjustment");
    FIELD_MAKE_OPTIONAL(calibrationStyle);
    FIELD (calibrationMaturity, "calibration Maturity (eg 10Y) for "
                  "quanto adjustment");
    FIELD_MAKE_OPTIONAL (calibrationMaturity);
}

CClassConstSP const ConvolutionEngine::TYPE = CClass::registerClassLoadMethod(
    "ConvolutionEngine", typeid(ConvolutionEngine), load);

void ConvolutionEngine::IIntoProduct::load(CClassSP& clazz){
    REGISTER_INTERFACE(IIntoProduct, clazz);
    EXTENDS(CModel::IModelIntoProduct);
}
    

CClassConstSP const ConvolutionEngine::IIntoProduct::TYPE = 
CClass::registerInterfaceLoadMethod(
    "ConvolutionEngine::IIntoProduct", typeid(IIntoProduct), load);

/** Returns par spread curves from a portfolio after being fetched
    from the market data cache */
class PortfolioSpreadCurvesAddin: public CObject {
public:
    // inner class for storing results for this addin - probably better done
    // with arrays rather than creating a whole new class just to store a
    // string and a DoubleArray
    class SpreadCurveOutput : public CObject {
    public:
        static CClassConstSP const TYPE;
        SpreadCurveOutput() :
            CObject(TYPE) {}
        SpreadCurveOutput(const string& nm, 
                          const DoubleArray& spds): 
            CObject(TYPE), name(nm), spreads(spds) {}
    private:
        static void load(CClassSP& clazz) {
            clazz->setPublic();
            REGISTER(SpreadCurveOutput, clazz);
            SUPERCLASS(CObject);
            EMPTY_SHELL_METHOD(defaultConstructor);

            //fields
            FIELD(name, "identifier");
            FIELD(spreads, "curve");
        }

        static IObject* defaultConstructor() {
            return new SpreadCurveOutput();
        }

        string      name;
        DoubleArray spreads;
    };

    typedef smartPtr<SpreadCurveOutput> SpreadCurveOutputSP;
    typedef array<SpreadCurveOutputSP, SpreadCurveOutput> SpreadCurveOutputArray;
    typedef smartPtr<SpreadCurveOutputArray> SpreadCurveOutputArraySP;

    static CClassConstSP const TYPE;

    IObjectSP run() {
        static string method = "PortfolioSpreadCurvesAddin::run";
        // clone instrument since we are altering it
        CInstrumentSP instCopy(inst.clone());
        // copy model as well since we might be filling it with market data
        ConvolutionEngineSP mdl(model.clone());
        //populate the model to get the credit index primarily
        mdl->getInstrumentAndModelMarket(market.get(), instCopy.get());
        // then create our product
        auto_ptr<ConvolutionProduct> prod(
            dynamic_cast<ConvolutionProduct*>(
               mdl->createProduct(instCopy.get())));

        if (!prod.get()) {
            throw ModelException(method, "Internal error: the product is not "
                                 "of type 'ConvolutionProduct' - This may be "
                                 "the case when this addin is requested on a "
                                 "Generalised CDO where the underlying loss "
                                 "config is not a tranche: only tranches "
                                 "support this addin at present.");
        }

        int numNames = prod->numNames();
        SpreadCurveOutputArraySP curves(new SpreadCurveOutputArray(numNames));

        //extract the curves
        for (int i = 0; i < numNames; i++) {
            SingleCreditAssetConstSP asset(prod->nameAsset(i));
            ICDSParSpreadsConstSP curve(asset->getParSpreadCurve());
            (*curves)[i].reset(new SpreadCurveOutput(
                                   curve->getName(),
                                   *curve->getParSpreads()));
        }
        // and return
        return curves;
    }
private:
    static void load(CClassSP& clazz) {
        REGISTER(PortfolioSpreadCurvesAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultConstructor);
        //fields
        FIELD(inst, "cdo instrument");
        FIELD(model, "model to fetch with");
        FIELD(market, "market data cache");

        Addin::registerObjectMethod(
            "PORTFOLIO_SPREAD_CURVES",
            Addin::RISK,
            "Returns par spread curves from a CDO "
            "portfolio after being fetched from the market data cache",
            false, //require a handle name
            Addin::expandSimple,
            &PortfolioSpreadCurvesAddin::run);
    }

    static IObject* defaultConstructor() {
        return new PortfolioSpreadCurvesAddin();
    }

    //for reflection
    PortfolioSpreadCurvesAddin(): CObject(TYPE) {}

    //fields
    CInstrumentSP        inst;   //the cdo instrument containing a portfolio
    ConvolutionEngineSP  model;  //the model
    CMarketDataSP        market; //the market data cache
};

// typedef works around VC71 bug

typedef PortfolioSpreadCurvesAddin::SpreadCurveOutput PortfolioSpreadCurvesAddin_SpreadCurveOutput;

CClassConstSP const PortfolioSpreadCurvesAddin_SpreadCurveOutput::TYPE = 
CClass::registerClassLoadMethod(
    "PortfolioSpreadCurvesAddin::SpreadCurveOutput",
    typeid(PortfolioSpreadCurvesAddin_SpreadCurveOutput), load);

// typedef works around VC71 bug
typedef PortfolioSpreadCurvesAddin::SpreadCurveOutputArray PortfolioSpreadCurvesAddin_SpreadCurveOutputArray;
DEFINE_TEMPLATE_TYPE_WITH_NAME("PortfolioSpreadCurvesAddin::SpreadCurveOutputArray", PortfolioSpreadCurvesAddin_SpreadCurveOutputArray);

CClassConstSP const PortfolioSpreadCurvesAddin::TYPE = 
CClass::registerClassLoadMethod(
    "PortfolioSpreadCurvesAddin", typeid(PortfolioSpreadCurvesAddin), load);


DRLIB_END_NAMESPACE

