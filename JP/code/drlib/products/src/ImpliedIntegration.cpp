//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ImpliedIntegration.cpp
//
//   Description : Implied Integration Algorithm 
//
//   Author      : Andrew J Swain
//
//   Date        :4 November 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ImpliedIntegration.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/MarketDataFetcherLNSpline.hpp"
#include "edginc/Control.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"

DRLIB_BEGIN_NAMESPACE

/** Constructor takes type of vol to use */
ImpliedIntegration::ImpliedIntegration(const string& volType):
    CModel(TYPE), midSteps(0), tailSteps(0), tailStdDev(0.0) {}


/** Override default createMDF in order to set the right MDF */
MarketDataFetcherSP ImpliedIntegration::createMDF() const{
    return MarketDataFetcherSP(new MarketDataFetcherLNSpline(volType));
}

IModel::WantsRiskMapping ImpliedIntegration::wantsRiskMapping() const {
    return riskMappingIrrelevant;
}

/** calculate single price and store result in results */
void ImpliedIntegration::Price(CInstrument*  instrument, 
                               CControl*     control, 
                               CResults*     results){
    static const string method = "ImpliedIntegration::Price";
    if (!IIntoProduct::TYPE->isInstance(instrument)){
        throw ModelException(method, "Instrument of type "+
                             instrument->getClass()->getName() +
                             " does not support ImpliedIntegration::IntoProduct");
    }
    IProduct* product = 0;
    try {
        if (instrument->priceDeadInstrument(control, results)) {
            return; // done for a dead instrument
        }
        
        // are we pricing or greeks ?
        isPricing = control ? control->isPricing() : false;

        // are we computing delta?
        if (!control)
        {
            isDeltaShift = false;
        }
        else
        {
            SensitivitySP sens = control->getCurrentSensitivity();
            if (!!sens && Delta::TYPE->isInstance(sens.get())) {
                isDeltaShift = true;
            } else {
                isDeltaShift = false;
            }
        }

        // cast to ImpliedIntegration::IIntoProduct
        IIntoProduct& intoProd = dynamic_cast<IIntoProduct&>(*instrument);
        // create the product
        product = intoProd.createProduct(this);

        product->price(this, control, results);
    } 
    catch (exception& e) {
        delete product;
        throw ModelException(e, method);
    }
    delete product;
}

/** integrate payoff */
double ImpliedIntegration::integrate(
    ImpliedIntegration::IProduct* product) const {
    static const string method = "ImpliedIntegration::integrate";
    try {
        PDFParams params(midSteps,
                         tailSteps,
                         PDFParams::default_numMiddleStrikes,
                         PDFParams::default_numTailStrikes,
                         1.0,
                         PDFParams::default_numMiddleStdDevs,
                         tailStdDev,
                         PDFParams::default_maxStdDevs,
                         PDFParams::default_maxFineStdDevs,
                         PDFParams::default_failProbs);

        DateTime        when = product->time();

        LinearImpliedSamplerSP sampler;

        SamplerMap::const_iterator iter = samplerMap.find(when);
        // 1) initial pricing => get LinearImpliedSampler from product and store it to reuse for greeks
        // 2) from outside (eg vdax) => same procedure
        // 3) create map for fwd starting products -- one sampler for "short maturity" and "long maturity"
        if (isPricing || iter == samplerMap.end()) {
            // create new sampler and append it to the map 
            sampler = LinearImpliedSamplerSP(product->sampler(when, params));
            samplerMap[when] = sampler;
        } else if (isTimeShift || isDeltaShift) {
            // also, if rolling through time, best to start from scratch
            // currently also do this for delta until same strikes approach is sorted out
            sampler = LinearImpliedSamplerSP(product->sampler(when, params));
        }
        else {
            // otherwise we want to keep the same strikes as last time
            // to avoid model noise
            LinearImpliedSamplerSP tweakSampler(product->sampler(when,params));

            sampler = LinearImpliedSamplerSP(new LinearImpliedSampler(tweakSampler.get(),
                                                                      (*iter).second->getStrikes(),
                                                                      (*iter).second->getPartition()));
        }

        LinearImpliedSampleArraySP sample(sampler->getImpliedSamples());

        Function1DDoubleSP payoff(product->payoff());
        Function1DDoubleSP integral(product->indefiniteIntegral());

        double i = (*sample)[0]->expectedValue(*payoff, *integral);

        return i;
    } 
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** what are the integration limits ? */
StrikesPartition* ImpliedIntegration::limits(
    ImpliedIntegration::IProduct* product) const {
    static const string method = "ImpliedIntegration::limits";
    try {
        PDFParams params(midSteps,
                         tailSteps,
                         PDFParams::default_numMiddleStrikes,
                         PDFParams::default_numTailStrikes,
                         1.0,
                         PDFParams::default_numMiddleStdDevs,
                         tailStdDev,
                         PDFParams::default_maxStdDevs,
                         PDFParams::default_maxFineStdDevs,
                         PDFParams::default_failProbs);

        LinearImpliedSamplerSP sampler(product->sampler(product->time(),params));

        StrikesPartitionSP partition(sampler->getPartition());

        return partition.release();
    } 
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** for VEGA_MATRIX - chose some semi-arbitrary strikes to report against */
/** mid/tailWidth -> determine hwo many strikes are used in mid/tail 
    regions -> get steps/width strikes
*/
DoubleArraySP ImpliedIntegration::sensitiveStrikes(
    StrikesPartition* bounds,
    int               midWidth,
    int               tailWidth) {
    static const string method = "ImpliedIntegration::sensitiveStrikes";
    try {
        int tailStrikes = tailSteps/tailWidth;  
        int midStrikes  = midSteps/midWidth;
        int tailStrikesDown = tailStrikes + 3;
 
        if (midStrikes < 2) {
            throw ModelException(method,
                                 "need at least 2 strikes in mid region");
        }
        if (tailStrikes < 2) {
            throw ModelException(method,
                                 "need at least 2 strikes in tail region");
        }
       
        double midHi  = bounds->midBrackets->upStrikes[0];
        double midLo  = bounds->midBrackets->downStrikes[0];

        int i;
        DoubleArray sensMid(midStrikes);
        for (i = 0; i < midStrikes; i++) {
            sensMid[i] = midLo + i *(midHi - midLo)/midStrikes;
        }

        // take the upper and lower strike spacings from the mid region to
        // determine the strike spacing in the upper and lower tails
        double upperSpacing = sensMid[midStrikes-1]-sensMid[midStrikes-2];
        DoubleArray sensUpper(tailStrikes);
        for (i = 0; i < tailStrikes; i++) {
            sensUpper[i]= midHi + i*upperSpacing;
        }

        // trickier at lower end as we're bounded at zero
        double lowerSpacing = sensMid[1]-sensMid[0];
        DoubleArray sensLower(0);
        for (i = 0; i < tailStrikesDown; i++) {
            double k = sensMid[0] - (tailStrikesDown-i)*lowerSpacing;
            if (Maths::isPositive(k)) {
                sensLower.push_back(k);
            }
        }
       
        // and stick them all together
        DoubleArraySP sensStrikes(new DoubleArray(0));
        for (i = 0; i < sensLower.size(); i++) {
            sensStrikes->push_back(sensLower[i]);
        }
        for (i = 0; i < sensMid.size(); i++) {
            sensStrikes->push_back(sensMid[i]);
        }
        for (i = 0; i < sensUpper.size(); i++) {
            sensStrikes->push_back(sensUpper[i]);
        }
         
        return sensStrikes;
    } 
    catch (exception& e) {
        throw ModelException(e, method);
    }
}
                                      

DoubleArraySP ImpliedIntegration::sensitiveStrikes(
    ImpliedIntegration::IProduct* product,
    int                           midWidth,
    int                           tailWidth) {
    static const string method = "ImpliedIntegration::sensitiveStrikes";
    try {
        // get the integration limits
        StrikesPartitionSP bounds(limits(product));

        return sensitiveStrikes(bounds.get(), midWidth, tailWidth);
    } 
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

DoubleArraySP ImpliedIntegration::sensitiveStrikes(const DateTime& today,
                                                   const DateTime& maturity,
                                                   const CAsset*   asset,
                                                   int             midWidth,
                                                   int             tailWidth) {
    static const string method = "ImpliedIntegration::sensitiveStrikes";
    try {
        PDFParams params(midSteps,
                         tailSteps,
                         PDFParams::default_numMiddleStrikes,
                         PDFParams::default_numTailStrikes,
                         1.0,
                         PDFParams::default_numMiddleStdDevs,
                         tailStdDev,
                         PDFParams::default_maxStdDevs,
                         PDFParams::default_maxFineStdDevs,
                         PDFParams::default_failProbs);

        LinearStrikeTSVolRequestSP volRequest(
            new LinearStrikeTSVolRequest(0.0,
                                         today,
                                         maturity,
                                         false)); 

        PDFRequestLNStrikeSP pdfRequest(new PDFRequestLNStrike(volRequest.get()));
        
        DateTimeArray datesTo(1, maturity);
        
        CAssetConstSP assetSP(asset);

        LinearImpliedSamplerSP sampler(new LinearImpliedSampler(assetSP,
                                                                volRequest,
                                                                today,
                                                                today,
                                                                datesTo,
                                                                params,
                                                                pdfRequest));

        StrikesPartitionSP bounds(sampler->getPartition());

        return sensitiveStrikes(bounds.get(), midWidth, tailWidth);
    } 
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

ImpliedIntegration::ImpliedIntegration():CModel(TYPE), isPricing(false),
                                         isTimeShift(false), isDeltaShift(false) {};

/** Does special things for Theta-type tweaks */
bool ImpliedIntegration::sensShift(Theta* shift) {
    // no changes for non-fwd starting things, however, what do we have to do for fwd starting?
    if (samplerMap.size()>0){
        // i.e. this is really is theta or similar, not the roll to now
        isTimeShift = true;
    }
    return true;
}

void ImpliedIntegration::sensRestore(Theta* shift) {
    isTimeShift = false;
}

/** implemented in order avoid problems when called several times for different dates */
void ImpliedIntegration::flush() {
    isPricing = true;
    isTimeShift = false;
    isDeltaShift = false;
    samplerMap.clear();
}

/** the field samplerMap is unregistered, hence, rewrite the clone method */
IObject* ImpliedIntegration:: clone() const {
    const string routine = "ImpliedIntegration::clone";        
    try {
        ImpliedIntegration* copy = dynamic_cast<ImpliedIntegration*>(CModel::clone());
        if(!copy) {
            throw ModelException(routine, "Clone method failed");
        }

        try {
            // Take deep copy of the map
            SamplerMap::const_iterator iter = samplerMap.begin();
            for(; iter != samplerMap.end(); ++iter) {
                LinearImpliedSamplerSP samplerClone(
                    dynamic_cast<LinearImpliedSampler*>((*iter).second->clone()));
                if(!samplerClone) {
                    throw ModelException("Failed to clone implied sampler for maturity " +
                                         (*iter).first.toString());
                }
                copy->samplerMap[(*iter).first] = samplerClone;
            }
        }
        catch (exception& ) {
            // don't die - give up on sampler map and proceed
            // greeks might be noisy but at least they don't fail
            samplerMap.clear();
        }

        return copy;
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


class ImpliedIntegrationHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ImpliedIntegration, clazz);
        SUPERCLASS(CModel);
        IMPLEMENTS(Theta::RestorableShift);
        EMPTY_SHELL_METHOD(defaultImpliedIntegration);
        FIELD(volType, "volType");
        FIELD(midSteps, "midSteps");
        FIELD(tailSteps, "tailSteps");
        FIELD(tailStdDev, "tailStdDev");

        FIELD(isPricing, "");
        FIELD_MAKE_TRANSIENT(isPricing);
        FIELD(isTimeShift, "");
        FIELD_MAKE_TRANSIENT(isTimeShift);
        FIELD(isDeltaShift, "");
        FIELD_MAKE_TRANSIENT(isDeltaShift);
}

    // for ImpliedIntegration::IIntoProduct
    static void loadIntoProduct(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER_INTERFACE(ImpliedIntegration::IIntoProduct, clazz);
        EXTENDS(Model::IModelIntoProduct);
    }

    static IObject* defaultImpliedIntegration(){
        return new ImpliedIntegration();
    }
};

CClassConstSP const ImpliedIntegration::TYPE = CClass::registerClassLoadMethod(
    "ImpliedIntegration", typeid(ImpliedIntegration), ImpliedIntegrationHelper::load);

bool  ImpliedIntegrationLoad() {
    return (ImpliedIntegration::TYPE != 0);
}

CClassConstSP const ImpliedIntegration::IIntoProduct::TYPE =
CClass::registerInterfaceLoadMethod("ImpliedIntegration::IIntoProduct",
                                    typeid(ImpliedIntegration::IIntoProduct), 
                                    ImpliedIntegrationHelper::loadIntoProduct);

DRLIB_END_NAMESPACE
