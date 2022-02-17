//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : CDOTrancheQuotes.cpp
//
//   Description : CDOTrancheQuotes is a market data container for CDO tranche quotes
//
//   Author      : Antoine Gregoire
//
//   Date        : December 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CDOTrancheQuotes.hpp"

DRLIB_BEGIN_NAMESPACE

/** Creates CDOTrancheQuotes with imput quotes (and empty market name) */
CDOTrancheQuotesSP CDOTrancheQuotes::create(
    const DateTime& valueDate,
    double lowStrike,
    double highStrike,
    ExpiryArraySP expiries,
    CDoubleArraySP spreads,
    CDoubleArraySP upfronts)
{
    // create new (empty) object
    CDOTrancheQuotesSP result(new CDOTrancheQuotes());
    
    // fill fields
    result->name = "";
    result->valueDate = valueDate;
    result->lowStrike = lowStrike;
    result->highStrike = highStrike;
    result->expiries = expiries;
    result->spreads = spreads;
    result->upfronts = upfronts;
    
    // validate
    result->validatePop2Object();
    
    // initialise maturityoExpiryIdxMap
    result->initMaturityoExpiryIdxMap();

    return result;
}

/** Destructor */
CDOTrancheQuotes::~CDOTrancheQuotes() {}

/** Returns the name of this object */
string CDOTrancheQuotes::getName() const {
	return name;
}

/** Populates the object with the market data that this object needs */
void CDOTrancheQuotes::getMarket(const IModel* model, const MarketData* market) {
    try {
        market->GetReferenceDate(valueDate);
        
        // initialise maturityoExpiryIdxMap
        initMaturityoExpiryIdxMap();
    } catch (exception& e){
        throw ModelException(e, "CDOTrancheQuotes::getMarket");
    }
}

/** Initiatlise maturityoExpiryIdxMap */
void CDOTrancheQuotes::initMaturityoExpiryIdxMap() {
    // clear map
    maturityoExpiryIdxMap.clear();    

    // build map
    for (int i = 0; i < expiries->size(); ++i) {
        maturityoExpiryIdxMap[(*expiries)[i]->toDate(valueDate)] = i;
    }
}

/** Called immediately after object constructed */
void CDOTrancheQuotes::validatePop2Object() {
    const string method = "CDOTrancheQuotes::validatePop2Object";
    try {
        // Check lowStrike < highStrike
        if (lowStrike >= highStrike) {
            throw ModelException(method, "Low strike greater than high strike.");
        }
        
        // Check expiries are in increasing order
        if (!Expiry::isIncreasing(valueDate, expiries.get())) {
            throw ModelException(method, "Expiries are not strictly increasing.");
        }
        
        // Check size of expiries and spreads arrays
        if (expiries->size() != spreads->size()) {
            throw ModelException(method, "Expiries and spreads arrays don't have the same size.");
        }

        // Check size of expiries and upfronts arrays (if any)
        if (upfronts.get() != 0) {
            if (expiries->size() != upfronts->size()) {
                throw ModelException(method, "Expiries and upfronts arrays don't have the same size.");
            }
        }
        
        // Check bid / ask data (if any)
        if (askSpreads.get() != 0 || bidSpreads.get() != 0) {
            // if ask (resp. bid) data is populated, we expect bid (resp. ask) data to be
            // populated as well
            if (askSpreads.get() == 0) {
                throw ModelException(method, "Bid spreads without ask spreads");
            }
            if (bidSpreads.get() == 0) {
                throw ModelException(method, "Ask spreads without bid spreads");
            }
            
            // check size of arrays
            int n = expiries->size();
            if (askSpreads->size() != n) {
                throw ModelException(method, "Expiries and ask spreads arrays don't have the same size.");
            }
            if (bidSpreads->size() != n) {
                throw ModelException(method, "Expiries and bid spreads arrays don't have the same size.");
            }
            if (askUpfronts.get() != 0) {
                if (expiries->size() != askUpfronts->size()) {
                    throw ModelException(method, "Expiries and ask upfronts arrays don't have the same size.");
                }
            }
            if (bidUpfronts.get() != 0) {
                if (expiries->size() != bidUpfronts->size()) {
                    throw ModelException(method, "Expiries and bid upfronts arrays don't have the same size.");
                }
            }
        }
    } catch (exception& e){
        throw ModelException(e, "CDOTrancheQuotes::validatePop2Object");
    }
}

/** Access to low strike */
double CDOTrancheQuotes::getLowStrike() const {
	return lowStrike;
}

/** Access to high strike */
double CDOTrancheQuotes::getHighStrike() const {
	return highStrike;
}

/** Access to expiries as DateTime */
DateTimeArraySP CDOTrancheQuotes::getMaturityDates() const {
    DateTimeArraySP result(new DateTimeArray(expiries->size()));
    for (int i = 0; i < expiries->size(); ++i) {
        (*result)[i] = (*expiries)[i]->toDate(valueDate);
    }
    return result;
}

/** Access to spreads */
double CDOTrancheQuotes::getSpread(const DateTime& maturity) const {
    return (*spreads)[getExpiryIdx(maturity)];
}

/** Access to upfronts */
double CDOTrancheQuotes::getUpfront(const DateTime& maturity) const {
    double upfront = 0.0;
    if (upfronts.get() != 0) {
        upfront = (*upfronts)[getExpiryIdx(maturity)];
    }
    return upfront;
}

/** Access to ask spreads */
double CDOTrancheQuotes::getAskSpread(const DateTime& maturity) const {
    if (askSpreads.get() != 0) {
        return (*askSpreads)[getExpiryIdx(maturity)];
    } else {
        throw ModelException("CDOTrancheQuotes::getAskSpread", "No ask spreads available");
    }
}

/** Access to ask upfronts */
double CDOTrancheQuotes::getAskUpfront(const DateTime& maturity) const {
    double askUpfront = 0.0;
    if (askUpfronts.get() != 0) {
        askUpfront = (*askUpfronts)[getExpiryIdx(maturity)];
    }
    return askUpfront;
}

/** Access to bid spreads */
double CDOTrancheQuotes::getBidSpread(const DateTime& maturity) const {
    if (bidSpreads.get() != 0) {
        return (*bidSpreads)[getExpiryIdx(maturity)];
    } else {
        throw ModelException("CDOTrancheQuotes::getBidSpread", "No bid spreads available");
    }
}

/** Access to bid upfronts */
double CDOTrancheQuotes::getBidUpfront(const DateTime& maturity) const {
    double bidUpfront = 0.0;
    if (bidUpfronts.get() != 0) {
        bidUpfront = (*bidUpfronts)[getExpiryIdx(maturity)];
    }
    return bidUpfront;
}

/** Returns true if bid / ask data (spreads and upfronts) is available */
bool CDOTrancheQuotes::hasBidAsk() const {
    return (askSpreads.get() != 0);
}

/**
 * Retrieve array index corresponding to "maturity",
 * throw an exception if not found
 * */
int CDOTrancheQuotes::getExpiryIdx(const DateTime& maturity) const {
    map<const DateTime, int>::const_iterator iter = maturityoExpiryIdxMap.find(maturity);
    if (iter != maturityoExpiryIdxMap.end()) {
        return iter->second;
    } else {
        throw ModelException("CDOTrancheQuotes::getExpiryIdx",
            "No quote available for maturity " + maturity.toString());
    }
}

/** Specific clone method to copy "maturityoExpiryIdxMap" field */
IObject* CDOTrancheQuotes::clone() const {
    CDOTrancheQuotes* copy = DYNAMIC_CAST(CDOTrancheQuotes, MarketObject::clone());
    copy->maturityoExpiryIdxMap = maturityoExpiryIdxMap;
    return copy;
}

/** Invoked when Class is 'loaded' */
void CDOTrancheQuotes::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(CDOTrancheQuotes, clazz);
    SUPERCLASS(MarketObject);
    EMPTY_SHELL_METHOD(defaultConstructor);
	
    FIELD(name, "Name");
    FIELD(lowStrike, "Attachment point");
    FIELD(highStrike, "Detachment point");
    FIELD(expiries, "Expiries");
    FIELD(spreads, "Mid spreads");
    FIELD(upfronts, "Mid upfront payments");
    FIELD_MAKE_OPTIONAL(upfronts);
    FIELD(bidSpreads, "Bid spreads");
    FIELD_MAKE_OPTIONAL(bidSpreads);
    FIELD(bidUpfronts, "Bid upfront payments");
    FIELD_MAKE_OPTIONAL(bidUpfronts);
    FIELD(askSpreads, "Ask spreads");
    FIELD_MAKE_OPTIONAL(askSpreads);
    FIELD(askUpfronts, "Ask upfront payments");
    FIELD_MAKE_OPTIONAL(askUpfronts);
    FIELD(valueDate, "Value date");
    FIELD_MAKE_TRANSIENT(valueDate);
}

/** Private constructor (only build instances of that class using reflection) */
CDOTrancheQuotes::CDOTrancheQuotes() :
	MarketObject(TYPE),
	upfronts(0), bidSpreads(0), bidUpfronts(0), askSpreads(0), askUpfronts(0) {}


/** Default constructor */
IObject* CDOTrancheQuotes::defaultConstructor() {
    return new CDOTrancheQuotes();
}

/** TYPE for CDOTrancheQuotes */
CClassConstSP const CDOTrancheQuotes::TYPE =
	CClass::registerClassLoadMethod(
    	"CDOTrancheQuotes",
    	typeid(CDOTrancheQuotes),
    	load);

/** TYPE for CDOTrancheQuotesArray */
DEFINE_TEMPLATE_TYPE(CDOTrancheQuotesArray);

/** TYPE for CDOTrancheQuotesWrapper */
DEFINE_TEMPLATE_TYPE(CDOTrancheQuotesWrapper);

/** TYPE for CDOTrancheQuotesWrapperArray */
DEFINE_TEMPLATE_TYPE(CDOTrancheQuotesWrapperArray);

DRLIB_END_NAMESPACE
