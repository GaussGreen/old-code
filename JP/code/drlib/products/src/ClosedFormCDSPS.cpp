//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : ClosedFormCDSPS.cpp
//
//   Description : Closed form model for taking CDS par spreads as input
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : August 30, 2001
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ClosedFormCDSPS.hpp"
#include "edginc/SRMFXDiffuse.hpp"
//#include "edginc/MarketDataFetcherCDS.hpp"
#include "edginc/MarketDataFetcherCIS.hpp"
#include "edginc/RhoPointwise.hpp"
#include "edginc/Control.hpp"
#include "edginc/IRGridPointCache.hpp"
#include "edginc/IInstrumentCollection.hpp"

DRLIB_BEGIN_NAMESPACE


ClosedFormCDSPS::ClosedFormCDSPS(): 
    Model(TYPE),
    // sensible default for ir vol interpolation (only used for
    // quanto adjustment)
    quantoCalibrationStyle("CMS"), 
    quantoCalibrationMaturity("10Y")
{}


ClosedFormCDSPS::ClosedFormCDSPS(CClassConstSP clazz):
    Model(clazz),
    // sensible default for ir vol interpolation (only used for
    // quanto adjustment)
    quantoCalibrationStyle("CMS"), 
    quantoCalibrationMaturity("10Y")
{}


/** Create a MarketDataFetcher which will be used for retrieving market data etc */
MarketDataFetcherSP ClosedFormCDSPS::createMDF() const {
    //return MarketDataFetcherSP(new MarketDataFetcherCDS(true));
    return MarketDataFetcherSP(new MarketDataFetcherCIS(true,false,false));//disable index basis
}

/** overridden to [reference] copy QuantoCDSAlgorithmSP */
IObject* ClosedFormCDSPS::clone() const{
    ClosedFormCDSPS* myCopy = DYNAMIC_CAST(ClosedFormCDSPS, Model::clone());
    myCopy->quantoAlgorithm = quantoAlgorithm; // just copy reference
    myCopy->irGridPtsCache = irGridPtsCache; // just copy reference
    return myCopy;
}

static void purifyFix(ClosedFormCDSPS::IProduct*     product,
                      ClosedFormCDSPS*               closedFormCDSPS,
                      CControl*                      control,
                      CResults*                      results){
    product->price(closedFormCDSPS, control, results);
}

/** calculate single price and store result in CResult */
void ClosedFormCDSPS::Price(CInstrument*  instrument, 
                            CControl*     control, 
                            CResults*     results){
    static const string method = "ClosedFormCDSPS::Price";
    IIntoProduct* intoProd;
    if (!IIntoProduct::TYPE->isInstance(instrument) ||
        !(intoProd = dynamic_cast<IIntoProduct*>(instrument))){
        throw ModelException(method, "Instrument of type "+
                             instrument->getClass()->getName() +
                             " does not support ClosedFormCDSPS::IntoProduct");
    }
    IProduct*     product = 0;
    try{
        if (control->isPricing() && quantoAlgorithm.get()){
            irGridPtsCache.reset(new IRGridPointCache(true));
            quantoAlgorithm->setIRGridPointCache(irGridPtsCache);
        }
        product = intoProd->createProduct(this);
        //product->price(this, control, results);
        purifyFix(product, this, control, results);
        if (control->isPricing() && quantoAlgorithm.get()){
            // switch off cache
            irGridPtsCache->setCachingMode(false);
        }
    } catch (exception& e){
        delete product;
        throw ModelException(e, method);
    }
    delete product;
}
    
/** Returns an instance of IAlgorithm. Typically the quantoCDSParSpreads
    parameter would be ignored but is there in case you want to switch
    the algorithm dependent upon some property of the quanto'd curve */
QuantoCDSParSpreads::IAlgorithmSP ClosedFormCDSPS::cdsQuantoAlgorithm(
    const QuantoCDSParSpreads* quantoCDSParSpreads) const{
    if (!quantoAlgorithm){
        createQuantoAlgorithm(false);
    }
    return quantoAlgorithm;
}

/** For quanto, need to adjust instEndDate to take into account IR spot vol
    calibration which looks at coupons of swap going beyond instEndDate */
DateTime ClosedFormCDSPS::endDate(const Sensitivity* sensitivity,
                                  const CInstrument* inst,
                                  const DateTime&    instEndDate) const {
    // we need a better way of ensuring that we cover any other similar
    // sensitivities to rho pointwise
    if (quantoAlgorithm.get() && RhoPointwise::TYPE->isInstance(sensitivity)){
        // only for quanto case
        OutputNameConstSP outputName =
            dynamic_cast<const RhoPointwise*>(sensitivity)->getMarketDataName();
        // juse reuse instEndDate as IR Vega pointwise end date
        return irGridPtsCache->rhoEndDate(outputName, instEndDate);
    }
    //return MaturityPeriod::toDate(30, "Y", instEndDate);
    return instEndDate;
}

/** Uses quanto calibration parameters to identify relevant points */
IRGridPointAbsArraySP ClosedFormCDSPS::getSensitiveIRVolPoints(
    OutputNameConstSP  outputName,
    const CInstrument* inst) const{
    if (!irGridPtsCache){
        throw ModelException("ClosedFormCDSPS::getSensitiveIRVolPoints",
                             "Internal error - no data for "+
                             outputName->toString());
    }
    return irGridPtsCache->sensitiveIRVolPoints(outputName);
}

IModel::WantsRiskMapping ClosedFormCDSPS::wantsRiskMapping() const {
    return riskMappingIrrelevant;
}

/** sets debug state to specified value. If true, then any calls to
    cdsQuantoAlgorithm will return an object that will cache debug data.
    This can be retrieved via getDebugInfo() */
void ClosedFormCDSPS::selectDebugState(bool switchOn){
    createQuantoAlgorithm(switchOn);
}

/** Returns debug info - may be null eg if the algorithm has not been 
    used */
IObjectSP ClosedFormCDSPS::getDebugInfo() const{
    if (!quantoAlgorithm){
        return IObjectSP();
    }
    return quantoAlgorithm->getDebugInfo();
}

/** Builds QuantoCDSAlgorithm object */
void ClosedFormCDSPS::createQuantoAlgorithm(bool debugOn) const{
    // decent error message in case empty strings have been passed in
    if (quantoCalibrationStyle.empty() || quantoCalibrationMaturity.empty()){
        throw ModelException("ClosedFormCDSPS:createQuantoAlgorithm", "Both "
                             "quantoCalibrationStyle and "
                             "quantoCalibrationMaturity model"
                             " parameters must be specified when cds curves "
                             "need to be quanto'd");
    }
    /* Note: the last 4 strings below are for the 'correlation swap'. The
       actual values seem to matter very little. The values here match
       maxim/aladdin.  */
    quantoAlgorithm.reset(new QuantoCDSAlgorithm("1F-STANDARD", // hard coded
                                                 quantoCalibrationStyle,
                                                 quantoCalibrationMaturity,
                                                 false, // skipIRBadVols?
                                                 SRMFXDiffuse::USE_LAST_LEVEL,
                                                 0.0, // fx cut off level
                                                 false, // floor negative vols?
                                                 // swap starting in 10Y
                                                 // maturity 10Y
                                                 "10Y", "10Y", "Act/365F", "A",
                                                 debugOn));
}

//------------------------------
// IHasForwardRatePricer methods
//------------------------------

/** Key method providing access to the pricer */
IForwardRatePricerSP ClosedFormCDSPS::getForwardRatePricer() const
{
    if (!forwardRateModel)
    {
        //not provided by the user, so supply the default
        return getDefaultForwardRatePricer();
    }
    else
    {
        //use the supplied pricer
        return forwardRateModel;
    }
}

//------------------------------
// CModel methods
//------------------------------
/** invoked after instruments has got its market data.  Allow model
 ** to get extra data */
void ClosedFormCDSPS::getMarket(const MarketData *market,
                                IInstrumentCollectionSP instruments)
{
    if(forwardRateModel.get())
    {
        IModel *model = this;
        forwardRateModel->getMarket(model, market);
    }

}



class ClosedFormCDSPSHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ClosedFormCDSPS, clazz);
        SUPERCLASS(Model);
        IMPLEMENTS(QuantoCDSParSpreads::IAlgorithmBuilder);
        EMPTY_SHELL_METHOD(defaultClosedFormCDSPS);
        FIELD(quantoCalibrationStyle,
                     "calibration style (eg CMS) for quanto"
                     " adjustment");
        FIELD_MAKE_OPTIONAL(quantoCalibrationStyle);
        FIELD(quantoCalibrationMaturity,
                     "calibration Maturity (eg 10Y) for "
                     "quanto adjustment");
        FIELD_MAKE_OPTIONAL(quantoCalibrationMaturity);
        FIELD       (forwardRateModel,    "A model capable of pricing all fees");
        FIELD_MAKE_OPTIONAL(forwardRateModel);
    }

    static IObject* defaultClosedFormCDSPS(){
        return new ClosedFormCDSPS();
    }
};

CClassConstSP const ClosedFormCDSPS::TYPE = CClass::registerClassLoadMethod(
    "ClosedFormCDSPS", typeid(ClosedFormCDSPS), ClosedFormCDSPSHelper::load);


bool ClosedFormCDSPSLoad() {
    return (ClosedFormCDSPS::TYPE != 0);
}


CClassConstSP const ClosedFormCDSPS::IIntoProduct::TYPE =
CClass::registerInterfaceLoadMethod("ClosedFormCDSPS::IIntoProduct",
                                    typeid(ClosedFormCDSPS::IIntoProduct), 0);


DRLIB_END_NAMESPACE
