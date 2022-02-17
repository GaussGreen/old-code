//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PseudoSimpleEquity.cpp
//
//   Description : pseudo equity
//
//   Date        : 30 Nov 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/PseudoSimpleEquity.hpp"
#include "edginc/VolRequestLN.hpp"
#include "edginc/Format.hpp"
#include "edginc/DividendCollector.hpp"


DRLIB_BEGIN_NAMESPACE

/** Calculate the settlement date associated with a given trade date */
DateTime PseudoSimpleEquity::settleDate(const DateTime& tradeDate) const{
    return parent->settleDate(tradeDate);
}

/** Returns the name (not the ISO code) of the yield curve used to
grow the stock */
string PseudoSimpleEquity::getYCName() const
{
    return parent->getYCName();
}

/** returns the spot price */
double PseudoSimpleEquity::getSpot() const{
    return parent->getEquity()->getPseudoSpot(horizonDate);
}

/** Returns fair value of stock price */
double PseudoSimpleEquity::fairValue() const
{
    // do I need this method ?
    throw ModelException("PseudoSimpleEquity::fairValue", "Not supported.");
    return parent->fairValue() / parent->getSpot() * getSpot();
}

/** returns the asset name */
string PseudoSimpleEquity::getName() const{
    return Format::toString("PSEUDO ") + parent->getName();
}

/** returns the real asset name */
string PseudoSimpleEquity::getTrueName() const{
    return parent->getName();
}

// the IMarketObservable interface for retrieving a single sample
double PseudoSimpleEquity::pastValue(const DateTime&       sampleDate,
                               const ObservationType*      obsType,
                               const ObservationSource*    source,
                               const FixingType*           fixType,
                               const IObservationOverride* overrides,
                               const SamplingConvention*   sampleRule) const{
    throw ModelException("PseudoSimpleEquity::pastValue", 
                "Cannot implement this method for PseudoSimpleEquity");
}

// IMarketObservable - retrieve a single observation date
// Returns false if obs is to be omitted
bool PseudoSimpleEquity::observationDate(const DateTime&           sampleDate,
                                         const ObservationSource*  source,
                                         const SamplingConvention* sampleRule,
                                         DateTime*                 obsDate) const {
    throw ModelException("PseudoSimpleEquity::observationDate", 
                "Cannot implement this method for PseudoSimpleEquity");
}

// the IMarketObservable interface for retrieving past samples events
double PseudoSimpleEquity::addPastSampleEvent(const DateTime&             sampleDate,
                                  const ObservationType*      obsType,
                                  const ObservationSource*    source,
                                  const FixingType*           fixType,
                                  const IObservationOverride* overrides,
                                  const SamplingConvention*   sampleRule,
                                  PastSamplesCollector*        collector) const {
    throw ModelException("PseudoSimpleEquity::addPastSampleEvent", 
                "Cannot implement this method for PseudoSimpleEquity");
}

// the IMarketObservable interface for 
// is the given date a holiday for the relevant source
bool PseudoSimpleEquity::isHoliday(const DateTime& sampleDate,
                             const ObservationSource*   source) const{
    throw ModelException("PseudoSimpleEquity::isHoliday", 
                "Cannot implement this method for PseudoSimpleEquity");
}

/** Calculates the expected spot price of the asset at the given date */
double PseudoSimpleEquity::fwdValue(const DateTime& date) const{
    CDoubleArray result(1);
    fwdValue(DateTimeArray(1, date),
             result);
    return result[0];
}

/** Calculates the expected spot prices of the asset at the given dates
    respecting any 'algorithmic' choices set in the FwdValueAlgorithm */
void PseudoSimpleEquity::fwdValue(const DateTimeArray&     dates,
                                  const FwdValueAlgorithm& algo,
                                  CDoubleArray&            result) const{
    if (algo.ignoreDivs()){
        throw ModelException("PseudoSimpleEquity::fwdValue", "ignoreDivs is not supported.");
    }
    fwdValue(dates,
             result);
}

/** Calculates the expected spot price of the asset at the given date if
    the spot price had the given value spot on spotDate */
double PseudoSimpleEquity::fwdFwd(const DateTime& spotDate,
                            double          spot, 
                            const DateTime& fwdDate) const {
    throw ModelException("PseudoSimpleEquity::fwdFwd", "Not supported.");
}

/** Calculates the stock floor at given dates assuming no dollar dividends 
    are paid passed horizonDate. In the absence of borrows and dividend yields 
    (continous and discrete) the stock floor is the PV of future dividends
    till horizonDate. */
void PseudoSimpleEquity::calcStockFloor(const DateTimeArray& dateList,
                                        CDoubleArray&        result) const {
    parent->getEquity()->calcStockFloor(dateList,
                                   horizonDate,
                                   result);
}

/** Calculates an array of forward prices which are assumed to be in 
    ascending order. The forward prices are those of the pseudo asset Z 
    that is associated with a dollar dividend paying equity S. See below for detail. */
void PseudoSimpleEquity::fwdValue(const DateTimeArray& dateList,
                                  CDoubleArray&        result) const {
    parent->getEquity()->calcPseudoFwdValue(dateList,
                                       horizonDate,
                                       result);
}

/** Given a vol base and a vol request, returns a processed vol where the 
    processed vol is that of the pseudo asset Z that is associated with a 
    dollar dividend paying equity S. The pseudo asset has the property that 
    Z = S - L where L is the floor of the stock price; it grows like S, 
    except for the dollar dividend drops. 
    NB Only simple equities and started options are supported at present. */
CVolProcessed* PseudoSimpleEquity::getProcessedVol(const CVolRequest* volRequest) const {
    static const string method = "PseudoSimpleEquity::getProcessedVol";
    try{
        const CVolRequestLN* volRequestLN = dynamic_cast<const CVolRequestLN*>(volRequest);
        if (!volRequestLN) {
            throw ModelException(method,
                                 "Only VolRequestLN types are supported; got "
                                 + volRequest->getClass()->getName() + ".");
        }
        return parent->getEquity()->getPseudoProcessedVol(parent->getVol().get(),
                                                     volRequestLN,
                                                     parent.get(),
                                                     horizonDate,
                                                     *divCritDates,
                                                     isCall,
                                                     noExerciseWindow);
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}                                                     

void PseudoSimpleEquity::getSensitiveStrikes(const CVolRequest* volRequest,
                                             OutputNameConstSP outputName,
                                             const SensitiveStrikeDescriptor& sensStrikeDesc,
                                             DoubleArraySP sensitiveStrikes) const {
    parent->getSensitiveStrikes(volRequest,
                                outputName,
                                sensStrikeDesc,
                                sensitiveStrikes);
}

/** The redirects everything to the SimpleEquity contained within */
bool PseudoSimpleEquity::accept(ICollector* collector) const{
    return parent->accept(collector);
}

/** constructor */
PseudoSimpleEquity::PseudoSimpleEquity(const EquityBase*  simpleEquity,
                                       const DateTimeArray& divCritDates,
                                       bool                 isCall,
                                       int                  noExerciseWindow):
CAsset(TYPE),
horizonDate(simpleEquity->getEquity()->getDivTransPeriodEndDate()),
divCritDates(copy(&divCritDates)),
isCall(isCall),
noExerciseWindow(noExerciseWindow){
    /* Convert the equity into an equity where the dollar dividends have been converted into
       dividend yields (where that applies. See Equity::convertDollarDivs for detail). */
    EquitySP eq(Equity::convertDollarDivs(*simpleEquity->getEquity()));
    parent = SimpleEquityConstSP(new SimpleEquity(eq.get(),
                                                  simpleEquity->getVol().get()));
}

PseudoSimpleEquity* PseudoSimpleEquity::create(const CAsset*        asset,
                                               const DateTimeArray& divCritDates,
                                               bool                 isCall,
                                               int                  noExerciseWindow){
    static const string method = "PseudoSimpleEquity::create";
    try{
        const IObject* obj = asset; // fix for gcc
        // only support simple equity
        const EquityBase* eq = dynamic_cast<const SimpleEquity*>(obj);
        if (!eq){ // try equityCache type
            eq = dynamic_cast<const EquityCache*>(obj);
        }
        if (!eq){
            throw ModelException(method, 
                                 "Need a SimpleEquity or EquityCache; got a "
                                 + asset->getClass()->getName() + ".");
        }
        return new PseudoSimpleEquity(eq,
                                      divCritDates,
                                      isCall,
                                      noExerciseWindow);
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

PDFCalculator* PseudoSimpleEquity::pdfCalculator(const PDFRequest* request) const {
    return parent->pdfCalculator(request);
}

/* for reflection */
PseudoSimpleEquity::PseudoSimpleEquity(): CAsset(TYPE){}
    
class PseudoSimpleEquityHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(PseudoSimpleEquity, clazz);
        SUPERCLASS(CAsset);
        EMPTY_SHELL_METHOD(defaultPseudoSimpleEquity);
        FIELD(horizonDate, "horizonDate");
        FIELD_MAKE_TRANSIENT(horizonDate);
        FIELD(divCritDates, "divCritDates");
        FIELD_MAKE_TRANSIENT(divCritDates);
        FIELD(parent, "SimpleEquity parent");
        FIELD_MAKE_TRANSIENT(parent);
        FIELD(isCall, "isCall");
        FIELD_MAKE_TRANSIENT(isCall);
        FIELD(noExerciseWindow, "noExerciseWindow");
        FIELD_MAKE_TRANSIENT(noExerciseWindow);
        ClassSetAcceptMethod(acceptDividendCollector);
    }

    /** pass on our dividends */
    static void acceptDividendCollector(const PseudoSimpleEquity* asset,
                                        DividendCollector*        collector){
        collector->processComponentAsset(asset->parent.get(), 1.0);
    }

    static IObject* defaultPseudoSimpleEquity(){
        return new PseudoSimpleEquity();
    }
};

CClassConstSP const PseudoSimpleEquity::TYPE = CClass::registerClassLoadMethod(
    "PseudoSimpleEquity", typeid(PseudoSimpleEquity), PseudoSimpleEquityHelper::load);

DRLIB_END_NAMESPACE

