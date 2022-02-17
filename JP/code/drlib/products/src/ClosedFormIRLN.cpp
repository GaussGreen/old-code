//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ClosedFormIRLN.cpp
//
//   Description : Closed Form Algorithm (asks instrument to do it)
//                 Interest rate flavour
//
//   Author      : Mark A Robson/Andrew J Swain
//
//   Date        : 28 February 2005
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ClosedFormIRLN.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/IRVol.hpp"
#include "edginc/Black.hpp"
#include "edginc/MarketDataFetcher.hpp"
#include "edginc/SensitiveIRVolPoints.hpp"
#include "edginc/MDFUtil.hpp"
#include "edginc/ModelFilter.hpp"

DRLIB_BEGIN_NAMESPACE

IObject* ClosedFormIRLN::clone() const {
    ClosedFormIRLN* copy = new ClosedFormIRLN();
    copy->smileStyle = smileStyle;
    return copy;
}

//// this class doesn't need to exist any more
class ClosedFormIRLN::MarketDataFetcherIR : public MarketDataFetcher {
public:
    ~MarketDataFetcherIR() {}

    MarketDataFetcherIR(const string &smileStyle): MarketDataFetcher(){
        // for requests for IRVolBase, get swaption vols (ie IRVol) - as
        // some products still have explicit IRVolBases in them
        setRetrievalMode(IRVolBase::TYPE, true, IRVolCommon::TYPE);
        if (smileStyle.empty()){
            MDFUtil::setUseSimpleIRVol(*this);
        } else {
            // get 2q smile data
            setRetrievalMode(IRCalib::SmileBase::TYPE, true,
                IRCalib::Smile2Q_TYPE());
            // get 1 factor ir parameters
            setRetrievalMode(IRCalib::Model::TYPE, true, 
                             IRCalib::getModelType(1));
            // get swaption vols when asked for ir vols if appropriate
            // This will let through IRCalibs unchanged (since they are not
            // derived from IRVol)
            setRetrievalMode(IRVolBase::TYPE, IYieldCurve::TYPE,
                             true, IRVolCommon::TYPE);
        }
    }
};


/**********************************************************************/
/*************** end of ClosedFormIRLN::MarketDataFetcherIR ***********/
/************************* begin of ClosedFormIRLN ********************/
/**********************************************************************/

/** calculate single price and store result in CResult */
void ClosedFormIRLN::Price(CInstrument*  instrument, 
                           CControl*     control, 
                           CResults*     results){
    if (!IIntoProduct::TYPE->isInstance(instrument)){
        throw ModelException("ClosedFormIRLN::Price", "Instrument of type "+
                             instrument->getClass()->getName() +
                             " does not support ClosedFormIRLN::IntoProduct");
    }
    IProduct*  product = 0;
    try{
        // cast to ClosedFormIRLN::IIntoProduct
        IIntoProduct& intoProd = dynamic_cast<IIntoProduct&>(*instrument);
        // create the product
        product = intoProd.createProduct(this);
        // and the invoke the pricing
        product->price(this, control, results);
    } catch (exception& e){
        delete product;
        throw ModelException(e, "ClosedFormIRLN::Price");
    }
    delete product;
}

/** Essentially relies on instrument implementing ISensitiveIRVolPoints.
    If not returns null. */
IRGridPointAbsArraySP ClosedFormIRLN::getSensitiveIRVolPoints(
    OutputNameConstSP outputName,
    const Instrument* inst) const{
    const ISensitiveIRVolPoints* vpImnt =
        dynamic_cast<const ISensitiveIRVolPoints*>(inst);
    return vpImnt? vpImnt->getSensitiveIRVolPoints(outputName, this): 
        IRGridPointAbsArraySP();
}

// a version of black that uses 2Q smile unless smileStyle
// is empty or set to NO_IR_Q_SMILE in which case it
// just does the (log)normal thing
double ClosedFormIRLN::black(
    bool             isCall, 
    double           fwd, 
    double           strike, 
    double           pv, 
    double           variance,
    const IRVolBase* vol) {
    static const string method("ClosedFormIRLN::black");
    try {
        if (smileStyle.empty() || smileStyle == NO_IR_Q_SMILE) {
            /* no smile */
            return Black::price(isCall, fwd, strike, pv, variance);
        }
        /* smile 2Q */
        IRCalib::SmileRequest smileRequest(smileStyle);
        if (!vol) {
            throw ModelException(method, "vol parameter is null");
        }
        CVolProcessedSP volProcessed(vol->getProcessedVol(&smileRequest,0));
        IRCalib::VolProcessed *volData = dynamic_cast<IRCalib::VolProcessed*>(volProcessed.get());
        if (!volData) {
            throw ModelException(method,"volProcessed should be of type IRCalib::VolProcessed");
        }

        const DoubleArray& smileParams = volData->getParams();

        if (smileParams.size()<3) {
            throw ModelException(method,"Number of ir vol smile params wrong");
        }

        // As per FI convention external meaning for QLeft/QRight was 0/0 for lognormal (baseline)
        // and 1/1 for normal. Internally these needed to be 1/1 and 0/0 respectively hence a 
        // rather unfortunate and confusing translation.
        return Black::price2Q(
            isCall, fwd, strike, pv, variance, 
            1.0 - smileParams[0], /* qLeft */
            1.0 - smileParams[1], /* qRight */
            smileParams[2] /* fwdShift */
        );
    } 
    catch (exception& e) { throw ModelException(e, method); }
}


/** Create a MarketDataFetcher which will be used for retrieving market data etc */
MarketDataFetcherSP ClosedFormIRLN::createMDF() const {
    MarketDataFetcherIR* mdf = new MarketDataFetcherIR(smileStyle);
    return MarketDataFetcherSP(mdf);
}

/** returns an ExposureHighlighter - a "model" that does everything except
actually price, so you get to see what market data it uses */
ExposureHighlighter* ClosedFormIRLN::exposureHighlighter() {
    return new IRVegaPointwiseExposureHighlighter(this);
}

IModel::WantsRiskMapping ClosedFormIRLN::wantsRiskMapping() const {
    return riskMappingIrrelevant;
}

ClosedFormIRLN::ClosedFormIRLN(): CModel(TYPE) {};

const string ClosedFormIRLN::NO_IR_Q_SMILE = "NO_IR_Q_SMILE";

class ClosedFormIRLNHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ClosedFormIRLN, clazz);
        SUPERCLASS(CModel);
        EMPTY_SHELL_METHOD(defaultClosedFormIRLN);
        //IMPLEMENTS(ISensitivePoints);
        FIELD(smileStyle, "smileStyle");
        FIELD_MAKE_OPTIONAL(smileStyle);
    }

    // for ClosedFormLN::IIntoProduct
    static void loadIntoProduct(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER_INTERFACE(ClosedFormIRLN::IIntoProduct, clazz);
        EXTENDS(Model::IModelIntoProduct);
    }

    static IObject* defaultClosedFormIRLN(){
        return new ClosedFormIRLN();
    }
};

CClassConstSP const ClosedFormIRLN::TYPE = CClass::registerClassLoadMethod(
    "ClosedFormIRLN", typeid(ClosedFormIRLN), ClosedFormIRLNHelper::load);
bool  ClosedFormIRLNLoad() {
    return (ClosedFormIRLN::TYPE != 0);
   }


CClassConstSP const ClosedFormIRLN::IIntoProduct::TYPE =
CClass::registerInterfaceLoadMethod("ClosedFormIRLN::IIntoProduct",
                                    typeid(ClosedFormIRLN::IIntoProduct), 
                                    ClosedFormIRLNHelper::loadIntoProduct);


DRLIB_END_NAMESPACE
