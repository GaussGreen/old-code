
#include "edginc/config.hpp"
#define QLIB_CORRELATION_CPP
#include "edginc/Correlation.hpp"
#include "edginc/Maths.hpp"
#include "edginc/RiskPropertySensitivity.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/Addin.hpp"
#include "edginc/FXAsset.hpp"
#include "edginc/MaturityTimePeriod.hpp"
#include "edginc/IModel.hpp"				// for CorrSwapBasis
#include "edginc/CorrelationCategory.hpp"	// for CorrSwapBasis
#include "edginc/TimeMetric.hpp"			// for CorrSwapBasis
#include "edginc/MarketData.hpp"			// for CorrSwapBasis
#include "edginc/Addin.hpp"
#include "edginc/MultiMarketFactors.hpp"
#include "edginc/NonPricingModel.hpp"
#include "edginc/MarketDataFetcherLNSpline.hpp"

DRLIB_BEGIN_NAMESPACE

const string Correlation::BENCHMARK_EXPIRY = "3Y";
const int Correlation::BENCHMARK_TIME = DateTime::END_OF_DAY_TIME;

/** Validation */
void Correlation::validatePop2Object(){
    static const string method("Correlation::validatePop2Object");
    if (name.empty()){
        name = asset1 < asset2? (asset1+ "_" + asset2): (asset2 + "_" + asset1);
    } else if (name == asset1 || name == asset2){
        throw ModelException(method, "Correlation's name must be different to"
                             " either asset's name");
    }
    if (fabs(correlation) > 1.0){
        throw ModelException( method, "Correlation "+name+
                              " must be between -1.0 and 1.0.");
    }
}

void Correlation::getMarket(const IModel* model, const MarketData* market){
    static const string method("Correlation::getMarket");
    try {
        /* need to fetch corrswap market data, if required -> ask MDF */
	    MarketDataFetcherSP mdf(model->getMDF());
	    DateTimeSP maturity(mdf->getCorrSwapExpiry());                 
        if (!(!maturity.get())) {
            /** set up a dummy metric */
		    string name = "Holiday";
		    const HolidaySP holiday(new Holiday(name, DateTimeArray(0), false));
		    TimeMetricSP metric(new TimeMetric(0.05, holiday.get()));
		    DateTime valueDate;
		    market->GetReferenceDate(valueDate);
		    timeDiff = metric->yearFrac(valueDate,*maturity.get());	

            string region1, region2; 
		    if (market->hasData(asset1, CorrelationCategory::TYPE)) { // if not availalbe -> basis = 0.0 (zeroObj)
                MarketObjectSP mo1cat = 
                    model->GetMarket(market, asset1, CorrelationCategory::TYPE);		
                CorrelationCategorySP category1 = CorrelationCategorySP::dynamicCast(mo1cat);
                if (!category1) {
                    throw ModelException(method, "Required correlation category object not retrievable");
                }
			    region1 = category1->getRegionName(); 
                if (region1.empty()) {
                    throw ModelException(method, "Entry regionName for CorrelationCategory "
                        + category1->getName() + " missing, but expected.");
                }
                MarketObjectSP mo1basisAdj = 
                    model->GetMarket(market, region1, CorrSwapBasisAdj::TYPE);
                basisAdj1 = CorrSwapBasisAdjSP::dynamicCast(mo1basisAdj);
                if (!basisAdj1) {
                    throw ModelException(method, "Required corr swap basis adj object not retrievable");
                }
                basisAdj1->getMarket(model,market);
            } else {
                basisAdj1 = CorrSwapBasisAdjSP(new CorrSwapBasisAdj(true));
            }

		    if (market->hasData(asset2, CorrelationCategory::TYPE)) { // if not available -> basis = 0.0 (zeroObj)
                MarketObjectSP mo2cat = 
                    model->GetMarket(market, asset2, CorrelationCategory::TYPE);
                CorrelationCategorySP category2 = CorrelationCategorySP::dynamicCast(mo2cat);
                if (!category2) {
                    throw ModelException(method, "Required correlation category object not retrievable");
                }
			    region2 = category2->getRegionName(); 
                if (region2.empty()) {
                    throw ModelException(method, "Entry regionName for CorrelationCategory "
                        + category2->getName() + " missing, but expected.");
                }
                MarketObjectSP mo2basisAdj = 
				    model->GetMarket(market, region2, CorrSwapBasisAdj::TYPE);
			    basisAdj2 = CorrSwapBasisAdjSP::dynamicCast(mo2basisAdj);
                if (!basisAdj2) {
                    throw ModelException(method, "Required corr swap basis adj object not retrievable");
                }
                basisAdj2->getMarket(model,market);
            } else {
                basisAdj2 = CorrSwapBasisAdjSP(new CorrSwapBasisAdj(true));
            }    

            if ( market->hasData(asset1, CorrelationCategory::TYPE)
                && market->hasData(asset2, CorrelationCategory::TYPE) 
                && (!CString::equalsIgnoreCase(region1, region2)) ) {            
			    const string& samplingAdjName = 
				    market->getCorrSwapSamplingAdjName(region1, region2);
			    MarketObjectSP moSamplingAdj = 
				    model->GetMarket(market, samplingAdjName, CorrSwapSamplingAdj::TYPE);
			    samplingAdj = CorrSwapSamplingAdjSP::dynamicCast((moSamplingAdj));
                if(!samplingAdj) {
                    throw ModelException(method, "Required corr swap sampling adj object not retrievable");
                }
                samplingAdj->getMarket(model, market);
            } else { 
                samplingAdj = CorrSwapSamplingAdjSP(new CorrSwapSamplingAdj(true));
            }
        } /*else { // create zero objects everywhere
            basisAdj1 = CorrSwapBasisAdjSP(new CorrSwapBasisAdj(true));
            basisAdj2 = CorrSwapBasisAdjSP(new CorrSwapBasisAdj(true));
            samplingAdj = CorrSwapSamplingAdjSP(new CorrSwapSamplingAdj(true));
        }*/
    } catch (exception& e){
        throw ModelException(e, method);
    }
}
       
double Correlation::getCorrelation() const{
    static const string method("Correlation::getCorrelation");
    try {
        /** zero objects should always be in place */
        if ( (!(basisAdj1.get())) || (!(basisAdj2.get())) || (!(samplingAdj.get())) ){
            throw ModelException(method, "Error in Correlation::getCorrelation");
        }
        double basisAdjustment1 = basisAdj1->getSqueezeForExpiry(timeDiff);
        double basisAdjustment2 = basisAdj2->getSqueezeForExpiry(timeDiff);
        double samplingAdjustment = samplingAdj->getSqueeze();    
        double adjCorrel  = correlation * (1.0 + 0.5 * (basisAdjustment1 + basisAdjustment2) + samplingAdjustment);
        return min(max(adjCorrel, -1.0),1.0); // COLLAR        
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

ExpirySP Correlation::getCorrExpiry() const{
    return corrExpiry;
}

// Sensitivity methods
/** Configure this correlation object so that under tweaking it behaves
    properly when it is a correlation between object of type clazz1 and
    an object of type clazz2 */
void Correlation::configureForSensitivities(CClassConstSP clazz1,
                                            CClassConstSP clazz2){
    sensType = NONE;
    // not perfect - but will do for the moment
    bool isCAsset1 = CAsset::TYPE->isAssignableFrom(clazz1);
    bool isCAsset2 = CAsset::TYPE->isAssignableFrom(clazz2);
    bool isFXAsset1 = isCAsset1 && FXAsset::TYPE->isAssignableFrom(clazz1);
    bool isFXAsset2 = isCAsset2 && FXAsset::TYPE->isAssignableFrom(clazz2);
    if (isCAsset1 && isCAsset2){
        if (isFXAsset1 || isFXAsset2){ // NB was originally exclusive or
            // this is suspect since FX v FX is apparently called FX_PHI(?)
            sensType = FX_PHI;
        }
        if (!isFXAsset1 && !isFXAsset2){
            sensType = PHI;
        }
    }
}

#define IAMSENSITIVE(CORREL) (sensType == PHI ? (CORREL)->asset :   \
                              sensType == FX_PHI ? (CORREL)->fx :   \
                              (CORREL)->other)

/** Returns true if this correlation is [really] sensitive to the
    supplied sensitivity */
bool Correlation::isSensitiveTo(const IPerNameSensitivity* sens) const{

    const RiskPropertySensitivity<Void>* rps =
        dynamic_cast<const RiskPropertySensitivity<Void> *>(sens);

    if (rps) {
        IScalarRiskPropertyConstSP p = rps->property();
        const RiskProperty<Correl>* prop =
            dynamic_cast<const RiskProperty<Correl>*>(p.get());
        if (prop) {
            ASSERT(!!prop->tag());
            return IAMSENSITIVE(prop->tag());
        }
    }

    return false;
}

/** Returns the name of the protected equity - used to determine
    whether to tweak the object */
string Correlation::sensName(const Correl* tag) const{
	ASSERT(tag);
    return IAMSENSITIVE(tag) ? getName() : "";
}

TweakOutcome Correlation::sensShift(const PropertyTweak<Correl>& tweak) {
    try {
        ASSERT(!!tweak.tag);
        if (tweak.tag->op == Correl::SQUEEZE) {
            correlation += (1 - correlation) * tweak.coefficient;
            return TweakOutcome(tweak.coefficient, false);
        }
        else {
            double d = (tweak.tag->op == Correl::ABSOLUTE ||
                           Maths::isNegative(correlation)) ? tweak.coefficient :
                                                             -tweak.coefficient;
            correlation += d;
            return TweakOutcome(d, false);
        }
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

string Correlation::getName() const{
    return name;
}

Correlation::Correlation(const string& name,
                         const string& nameAsset1,
                         const string& nameAsset2,
                         double        correlation):
    CorrelationCommon(Correlation::TYPE, name, nameAsset1, nameAsset2), 
    correlation(correlation), 
    corrExpiry(ExpirySP(new MaturityTimePeriod(BENCHMARK_EXPIRY, BENCHMARK_TIME))), 
    basisAdj1(CorrSwapBasisAdjSP(new CorrSwapBasisAdj(true))),
    basisAdj2(CorrSwapBasisAdjSP(new CorrSwapBasisAdj(true))),
    timeDiff(0.0),
    samplingAdj(CorrSwapSamplingAdjSP(new CorrSwapSamplingAdj(true))),
    sensType(NONE) {
    validatePop2Object();
}

Correlation::Correlation():CorrelationCommon(TYPE), 
                           correlation(0.0), 
                           corrExpiry(ExpirySP(new MaturityTimePeriod(BENCHMARK_EXPIRY, BENCHMARK_TIME))), 
                           basisAdj1(CorrSwapBasisAdjSP(new CorrSwapBasisAdj(true))),
                           basisAdj2(CorrSwapBasisAdjSP(new CorrSwapBasisAdj(true))),
                           timeDiff(0.0),
                           samplingAdj(CorrSwapSamplingAdjSP(new CorrSwapSamplingAdj(true))),
                           sensType(NONE){}

/** Invoked when Class is 'loaded' */
void Correlation::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(Correlation, clazz);
    SUPERCLASS(CorrelationCommon);
    IMPLEMENTS(ITweakableWithRespectTo<Correl>);
    EMPTY_SHELL_METHOD(defaultCorrelation);
    clazz->enableCloneOptimisations();
    FIELD(name, "name for correlation");
    FIELD_MAKE_OPTIONAL(name);
    FIELD(asset1,      "Name of Asset 1");
    FIELD(asset2,      "Name of Asset 2");
    FIELD(correlation, "Correlation");
    FIELD(corrExpiry, "ref date for corr");
    FIELD_MAKE_OPTIONAL(corrExpiry);
    FIELD_NO_DESC(sensType);
    FIELD_MAKE_TRANSIENT(sensType);
	FIELD(samplingAdj, "samplingAdj");
	FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(samplingAdj);
    FIELD(basisAdj1, "basisAdj1");
	FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(basisAdj1);
    FIELD(basisAdj2, "basisAdj2");
	FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(basisAdj2);
    FIELD(timeDiff, "for corr swap maturity");
    FIELD_MAKE_TRANSIENT(timeDiff);
}

IObject* Correlation::defaultCorrelation(){
    return new Correlation();
}

CClassConstSP const Correlation::TYPE = CClass::registerClassLoadMethod(
    "Correlation", typeid(Correlation), load);

// initialise type for array of Correlations
DEFINE_TEMPLATE_TYPE(CorrelationArray);

class AddCorrelationAddin: public CObject{
    static CClassConstSP const TYPE;

    /* addin parameters - MarketData cache & data to add*/
    CMarketDataSP     market;
    StringArray       vertAxis;
    StringArray       horAxis;
    DoubleMatrix      correlations;

    /** clone market data cache and add to it */
    static IObjectSP add(AddCorrelationAddin *params){
        static const string routine("AddCorrelationAddin::add");
        CMarketDataSP market(copy(params->market.get()));

        if (params->vertAxis.size() != params->correlations.numRows() ||
            params->horAxis.size() != params->correlations.numCols()){
            throw ModelException(routine, "Mismatch between size of matrix and"
                                 " length of axes");
        }
        DoubleMatrix& matrix = params->correlations;
        map<string, double> corrNames; // keep track of what we've added
        for (int i = 0; i < matrix.numRows(); i++){
            for (int j = 0; j < matrix.numCols(); j++){
                MarketObjectSP object(new Correlation(
                    "", params->vertAxis[i], params->horAxis[j],
                    matrix[j][i]));
                // see if we've done this before
                const string& name = params->vertAxis[i] < params->horAxis[j]?
                    params->vertAxis[i]+"_"+params->horAxis[j]:
                    params->horAxis[j]+"_"+params->vertAxis[i];
                map<string, double>::const_iterator iter = 
                    corrNames.find(name);
                if (iter == corrNames.end()){
                    market->AddData(object);
                    corrNames[name] = matrix[j][i];
                } else {
                    // check values are the same
                    if (!Maths::equals(iter->second, matrix[j][i])){
                        throw ModelException(routine, "Inconsistent values for"
                                             " correlation between "+
                                             params->vertAxis[i] + " and "+
                                             params->horAxis[j]);
                    }
                }
            }
        }
        return market;
    }
 
    AddCorrelationAddin(): CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(AddCorrelationAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultAddCorrelationAddin);
        FIELD(market, "market data cache");
        FIELD(vertAxis, "Names for Vertical Axis");
        FIELD(horAxis, "Names for Horizontal Axis");
        FIELD(correlations, "The correlations");
        // register addin for adding objects to cache
        Addin::registerClassObjectMethod(
            "ADD_CORRELATIONS_TO_MARKET",
            Addin::MARKET,
            "copy market data and add correlation matrix to it",
            TYPE,
            true,
            Addin::returnHandle,
            (Addin::ObjMethod*)add);

    }

    static IObject* defaultAddCorrelationAddin(){
        return new AddCorrelationAddin();
    }
};

CClassConstSP const AddCorrelationAddin::TYPE= CClass::registerClassLoadMethod(
    "AddCorrelationAddin", typeid(AddCorrelationAddin), load);

DRLIB_END_NAMESPACE
