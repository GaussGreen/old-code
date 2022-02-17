//----------------------------------------------------------------------------
//
//   Group       : QR Equities London
//
//   Filename    : CorrSwapSamplingAdj.cpp
//
//   Description : 
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_CorrSwapSamplingAdj_CPP
#include "edginc/Atomic.hpp"
#include "edginc/CorrSwapSamplingAdj.hpp"
#include "edginc/Correlation.hpp"
#include "edginc/MarketDataFetcher.hpp"
#include "edginc/MarketData.hpp"


DRLIB_BEGIN_NAMESPACE


CorrSwapSamplingAdj::~CorrSwapSamplingAdj() {}

/** Constructor */
CorrSwapSamplingAdj::CorrSwapSamplingAdj() : MarketObject(TYPE), isZeroObj(false) {
    // none of the fields is optional,  the "name" field is built up in validatePop2Object ... 
}

/** Constructor for zero object */
CorrSwapSamplingAdj::CorrSwapSamplingAdj(bool isZeroObj) : MarketObject(TYPE), 
name(""), region1(""), region2(""), squeeze(0.0), isZeroObj(isZeroObj) {
    static const string method("CorrSwapSamplingAdj::CorrSwapSamplingAdj"); 
    try {
        if (!isZeroObj) {
            throw ModelException(method, "This constructor shall only be used to create a zero object.");
        }
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Validation */
void CorrSwapSamplingAdj::validatePop2Object() {
	static const string method("CorrSwapSamplingAdj::validatePop2Object");
	try {
        if (name.empty()) {
            name = (region1 < region2) ? (region1 + "_" + region2) : (region2 + "_" + region1);
        } else if ( CString::equalsIgnoreCase(name, region1) || CString::equalsIgnoreCase(name, region2) ) {
                throw ModelException(method, 
                    "Correl Swap Samlping Adj's name must be different toeither region's name");
        }    
        if (CString::equalsIgnoreCase(region1, region2)) {
            throw ModelException(method, "Corr swap sampling adj can only be defined cross regional, "
                "but is defined inter regional for " + region1);
        }
        if ( Maths::isPositive(fabs(squeeze) - 1.0) ) {
            throw ModelException(method, "Squeeze for corr swap sampling adjustment is "
                + Format::toString(squeeze)
                + " for "
                + name
                + " but must be between -1.0 and 1.0");
        }
    } catch (exception& e){
        throw ModelException(e, method);
    }
}
        
/** Returns the squeeze */ 
double CorrSwapSamplingAdj::getSqueeze() const {
	return squeeze;
}

/** Returns the CorrSwapSamplingAdj's name */
string CorrSwapSamplingAdj::getName() const {
	return name;
}

/** Tweak (sensitivity) support */
string CorrSwapSamplingAdj::sensName(const CorrSwapSamplingAdjAbsolute*) const {
    return name;
}

TweakOutcome CorrSwapSamplingAdj::sensShift(const PropertyTweak<CorrSwapSamplingAdjAbsolute>& tweak) {
    static const string method = "CorrSwapSamplingAdj::sensShift";
    try {
        if (isZeroObj) {
            return TweakOutcome(tweak.coefficient, false); 
        }
        if (!Maths::isZero(tweak.coefficient)) {
            squeeze += tweak.coefficient;            
        }
        return TweakOutcome(tweak.coefficient, true); 
    } catch (exception &e) {
        throw ModelException(e, method, "CorrSwapSamplingAdj tweak failed for " + getName());
    }
}

/** Initialises this piece of market data - 
records the pair of names idenitfying the correlation */
void CorrSwapSamplingAdj::initialise(MarketData* market) {
	string marketCorrSwapSamplingAdjName;
	bool   hasCorrSwapSamplingAdj = market->hasCorrSwapSamplingAdjData(region1, region2);

	if ( hasCorrSwapSamplingAdj ) {
		marketCorrSwapSamplingAdjName = market->getCorrSwapSamplingAdjName(region1, region2);
	}

	if ( !hasCorrSwapSamplingAdj || marketCorrSwapSamplingAdjName == name ) {
		market->setCorrSwapSamplingAdjName(region1, region2, name);
	} else {
		throw ModelException("CorrelationTerm::initialise", 
			"There is already a corrSwapSamplingAdj in the cache for assets " 
			+ region1 + " and " + region2 + " with name " + marketCorrSwapSamplingAdjName + ".\n"
			" Cannot add a corrSwapSamplingAdj for the same regions with a different name (" + name + ").");
	}
}


class CorrSwapSamplingAdjHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CorrSwapSamplingAdj, clazz);
        SUPERCLASS(MarketObject);
        IMPLEMENTS(ITweakableWithRespectTo<CorrSwapSamplingAdjAbsolute>);
        EMPTY_SHELL_METHOD(defaultCorrSwapSamplingAdj);
        FIELD(name, "name for correl swap basis");
        FIELD_MAKE_OPTIONAL(name);
        FIELD(region1, "name for region1"); 
        FIELD(region2, "name for region2"); 
        FIELD(squeeze, "squeeze");		
        FIELD(isZeroObj, "is zero obj");
        FIELD_MAKE_TRANSIENT(isZeroObj);
        // by default don't get CorrSwapSamplingAdj
        MarketDataFetcher::setDefaultRetrievalMode(CorrSwapSamplingAdj::TYPE,
                                                   false, NULL);
    }

    static IObject* defaultCorrSwapSamplingAdj(){
        return new CorrSwapSamplingAdj();
    }
};

CClassConstSP const CorrSwapSamplingAdj::TYPE = CClass::registerClassLoadMethod(
    "CorrSwapSamplingAdj", typeid(CorrSwapSamplingAdj), CorrSwapSamplingAdjHelper::load);

bool CorrSwapSamplingAdjLoad() {
    return (CorrSwapSamplingAdj::TYPE != 0);
}

DRLIB_END_NAMESPACE
