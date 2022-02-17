//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : CDOQuotes.cpp
//
//   Description : CDOQuotes is a market data container for index tranche quotes
//
//   Author      : Antoine Gregoire
//
//   Date        : December 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CDOQuotes.hpp"
#include "edginc/Format.hpp"
#include "edginc/TimePoint2D.hpp"
#include "edginc/BenchmarkDate.hpp"
#include <set>

DRLIB_BEGIN_NAMESPACE

/** TYPE for CDOQuotesWrapper */
DEFINE_TEMPLATE_TYPE(CDOQuotesWrapper);

/**
 * Generate a fee leg corresponding to the given
 * (lowStrike, highStrike, maturity) point.
 * Will fail if that point is not quoted.
 * */
ICreditFeeLegSP CDOQuotes::generateFeeLeg(
    double lowStrike, double highStrike, const DateTime& maturity) const
{
    // Retrieve tranche quotes index
    map<pair<double, double>, int>::const_iterator iter =
        strikesToTrancheQuotesIdxMap.find(pair<double, double>(lowStrike, highStrike));
    if (iter != strikesToTrancheQuotesIdxMap.end()) {
        int trancheQuoteIdx = iter->second;
        
        double spread = (*(*trancheQuotes)[trancheQuoteIdx])->getSpread(maturity);
        double upfront = (*(*trancheQuotes)[trancheQuoteIdx])->getUpfront(maturity);
        
        return trancheConventions->generateFeeLeg(
            trancheConventions->startDate(valueDate),
            maturity,
            spread,
            upfront);
    } else {
        throw ModelException(
            "CDOQuotes::generateFeeLeg",
            "No quotes available for tranche " +
            Format::toString(lowStrike) +
            "-" +
            Format::toString(highStrike));
    }
}

/**
 * Generate a fee leg using this CDOQuotes convention but with
 * overriden spread and upfront payment
 * */
ICreditFeeLegSP CDOQuotes::generateFeeLegOverride(
    double spread, double upfront, const DateTime& maturity) const
{
    return trancheConventions->generateFeeLeg(
        trancheConventions->startDate(valueDate),
        maturity,
        spread,
        upfront);
}

/**
 * Generate a contingent leg corresponding to the given
 * (lowStrike, highStrike, maturity) point.
 * Will fail if that point is not quoted.
 * */
ICreditContingentLegSP CDOQuotes::generateContingentLeg(
    double lowStrike, double highStrike, const DateTime& maturity) const
{
    // Retrieve tranche quotes index
    map<pair<double, double>, int>::const_iterator iter =
        strikesToTrancheQuotesIdxMap.find(pair<double, double>(lowStrike, highStrike));
    if (iter != strikesToTrancheQuotesIdxMap.end()) {
        int trancheQuoteIdx = iter->second;
        
        double spread = (*(*trancheQuotes)[trancheQuoteIdx])->getSpread(maturity);
        double upfront = (*(*trancheQuotes)[trancheQuoteIdx])->getUpfront(maturity);
        
        return trancheConventions->generateContingentLeg(
            trancheConventions->startDate(valueDate),
            maturity,
            spread,
            upfront);
    } else {
        throw ModelException(
            "CDOQuotes::generateContingentLeg",
            "No quotes available for tranche " +
            Format::toString(lowStrike) +
            "-" +
            Format::toString(highStrike));
    }
}

/** access to recover notional flag in creditLegConventions */
bool CDOQuotes::getRecoverNotional() const
{
    return trancheConventions->getRecoverNotional();
}

/**
 * Generate a contingent leg using this CDOQuotes convention but with
 * overriden spread and upfront payment
 * */
ICreditContingentLegSP CDOQuotes::generateContingentLegOverride(
    double spread, double upfront, const DateTime& maturity) const
{
    return trancheConventions->generateContingentLeg(
        trancheConventions->startDate(valueDate),
        maturity,
        spread,
        upfront);
}

/** Destructor */
CDOQuotes::~CDOQuotes() {}

/** Returns the name of this object */
string CDOQuotes::getName() const {
	return name;
}

/** Populates the object with the market data that this object needs */
void CDOQuotes::getMarket(const IModel* model, const MarketData* market) {
    try {
        market->GetReferenceDate(valueDate);

        // populates trancheQuotes
		for (int i = 0; i < trancheQuotes->size(); ++i) {
            (*(*trancheQuotes)[i]).getData(model, market);
		}

        trancheConventions->getMarket(model, market);
        
        // builds strikesToTrancheQuotesIdxMap
        initStrikesToTrancheQuotesIdxMap();
    } catch (exception& e){
        throw ModelException(e, "CDOQuotes::getMarket");
    }
}

/** Init strikesToTrancheQuotesIdxMap */
void CDOQuotes::initStrikesToTrancheQuotesIdxMap() {
    // clear map
    strikesToTrancheQuotesIdxMap.clear();

    // build strikesToTrancheQuotesIdxMap
    for (int i = 0; i < trancheQuotes->size(); ++i) {
        strikesToTrancheQuotesIdxMap[pair<double, double>(
            (*(*trancheQuotes)[i])->getLowStrike(),
            (*(*trancheQuotes)[i])->getHighStrike())] = i;
    }
}

/** Called immediately after object constructed */
void CDOQuotes::validatePop2Object() {
    try {
    	//TODO
    } catch (exception& e){
        throw ModelException(e, "CDOQuotes::validatePop2Object");
    }
}

/** Access to tranche quotes */
CDOTrancheQuotesArrayConstSP CDOQuotes::getTrancheQuotes() const {
	CDOTrancheQuotesArraySP result =
		CDOTrancheQuotesArraySP(new CDOTrancheQuotesArray(trancheQuotes->size()));
	for (int i = 0; i < trancheQuotes->size(); ++i) {
		(*result)[i] = (*trancheQuotes)[i]->getSP();
	}
	return result;
}

/** Access to portfolio */
CDOPortfolioConstSP CDOQuotes::getPortfolio() const {
    return portfolio;
} 

/** Access to discount curve name */
string CDOQuotes::getDiscountName() const {
    return trancheConventions->getDiscountName();
}

/** Access to discount curve */
YieldCurveConstSP CDOQuotes::getDiscount() const {
    return trancheConventions->getDiscount();
}

/** Access to value date */
const DateTime& CDOQuotes::getValueDate() const {
    return valueDate;
}

/** Specific clone method to copy "strikesToTrancheQuotesIdxMap" field */
IObject* CDOQuotes::clone() const {
    CDOQuotes* copy = DYNAMIC_CAST(CDOQuotes, MarketObject::clone());
    copy->strikesToTrancheQuotesIdxMap = strikesToTrancheQuotesIdxMap;
    return copy;
}

/**
 * Creates a new CDOQuotes object using the same conventions but
 * different tranche quotes.
 * The new CDOQuotes will have no "market name".
 * */
CDOQuotesSP CDOQuotes::createWithNewQuotes(
    ExpiryArraySP expiries,
    CDoubleArraySP lowStrikes,
    CDoubleArraySP highStrikes,
    CDoubleArraySP spreads,
    CDoubleArraySP upfronts) const
{
    static const string method = "CDOQuotes::createWithNewQuotes";
    int i, j;
    
    int n = lowStrikes->size();

    // initialise upfronts with 0.0 if not defined
    if (upfronts.get() == 0) {
        upfronts.reset(new CDoubleArray(n, 0.0));
    }

    // checks all array have the same length n
    if (highStrikes->size() != n || spreads->size() != n || upfronts->size() != n) {
        throw ModelException(method,
            "low strikes, high strikes, spreads and upfronts arrays do not have the same length"); 
    }
    
    // creates new CDOQuotes object with same conventions
    CDOQuotesSP newQuotes(DYNAMIC_CAST(CDOQuotes, MarketObject::clone()));
    
    // removes market name
    newQuotes->name = "";

    // defines a map from (low strike, high strike) to a maturity ordered set of (maturity, spread, upfront) points
    // so each entry of this map corresponds to a kind of "pseudo tranche quotes" that will be
    // converted into a "real" tranche quotes in a second step below
    typedef map<pair<double, double>, set<TimePoint2D, TimePoint2D::CompareDateC1C2> > StrikesToQuotesMap;
    StrikesToQuotesMap strikesToQuotesMap;
    
    // populates strikesToQuotesMap
    for (i = 0; i < n; ++i) {
		strikesToQuotesMap[pair<double, double>((*lowStrikes)[i], (*highStrikes)[i])].insert(
            TimePoint2D((*expiries)[i]->toDate(valueDate), (*spreads)[i], (*upfronts)[i]));
	}

    // transforms strikesToQuotesMap into CDOTrancheQuotesArray ("pseudo" -> "real" tranche quotes)
    
    int nbTrancheQuotes = strikesToQuotesMap.size();
    CDOTrancheQuotesArraySP newTrancheQuotes(new CDOTrancheQuotesArray(nbTrancheQuotes));

    double tqLowStrike, tqHighStrike; // "tq" stands for tranche quote
    ExpiryArraySP tqExpiries;
    CDoubleArraySP tqSpreads, tqUpfronts;
    int nbTQExpiries;
    
    // loops through "pseudo" tranche quotes
    StrikesToQuotesMap::iterator mapIter = strikesToQuotesMap.begin();
    for (i = 0; mapIter != strikesToQuotesMap.end(); mapIter++){
        tqLowStrike = (*mapIter).first.first;
        tqHighStrike = (*mapIter).first.second;
        nbTQExpiries = (*mapIter).second.size();
        tqExpiries.reset(new ExpiryArray(nbTQExpiries));
        tqSpreads.reset(new CDoubleArray(nbTQExpiries));
        tqUpfronts.reset(new CDoubleArray(nbTQExpiries));
        
        // loops through this "pseudo" tranche quote to populate expiries, spreads and upfronts
        set<TimePoint2D, TimePoint2D::CompareDateC1C2>::iterator setIter = (*mapIter).second.begin();
        for (j = 0; setIter != (*mapIter).second.end(); setIter++) {
			(*tqExpiries)[j].reset(new BenchmarkDate(setIter->getDate()));
            (*tqSpreads)[j] = setIter->getCoord1(); 
            (*tqUpfronts)[j] = setIter->getCoord2(); 

            j++;
		}
        
        // builds "real" tranche quote
        (*newTrancheQuotes)[i] = CDOTrancheQuotes::create(
            valueDate, tqLowStrike, tqHighStrike, tqExpiries, tqSpreads, tqUpfronts);

        i++;
	}
    
    // puts new tranche quotes into newQuotes
    newQuotes->trancheQuotes.reset(new CDOTrancheQuotesWrapperArray(newTrancheQuotes->size()));
    for (i = 0; i < newTrancheQuotes->size(); ++i) {
        // creates wrapper with empty name
        (*newQuotes->trancheQuotes)[i] = CDOTrancheQuotesWrapperSP(new CDOTrancheQuotesWrapper(""));
        
        // populates wrapper
        (*(*newQuotes->trancheQuotes)[i]).setObject((*newTrancheQuotes)[i]);
    }

    // builds strikesToTrancheQuotesIdxMap
    newQuotes->initStrikesToTrancheQuotesIdxMap();

    return newQuotes;
}

/** Invoked when Class is 'loaded' */
void CDOQuotes::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(CDOQuotes, clazz);
    SUPERCLASS(MarketObject);
    IMPLEMENTS(ICDOQuotesGenerator);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(name, "Name");
    FIELD(trancheQuotes, "Tranche quotes");
    FIELD(portfolio,
        "Portfolio (including single names betas and par spreads)");
    FIELD(trancheConventions, "Tranche conventions");
    FIELD(valueDate, "Value date");
    FIELD_MAKE_TRANSIENT(valueDate);
}

/** Private constructor (only build instances of that class using reflection) */
CDOQuotes::CDOQuotes() : MarketObject(TYPE), portfolio(0), trancheConventions(0) {}

/**
 * Public Constructor.
 * Creates CDOQuotes object explicitely from
 * a flattened trancheQuotes
 * */
CDOQuotes::CDOQuotes(
    DateTime       valueDate,
    CDOPortfolioSP   portfolio,
    ICreditLegConventionSP trancheConventions,
    ExpiryArraySP  expiries,
    CDoubleArraySP lowStrikes,
    CDoubleArraySP highStrikes,
    CDoubleArraySP spreads,
    CDoubleArraySP upfronts)
: MarketObject(TYPE), portfolio(portfolio), trancheConventions(trancheConventions), valueDate(valueDate)
{
    static const string method = "CDOQuotes::CDOQuotes";
    int i, j;
    
    int n = lowStrikes->size();

    // initialise upfronts with 0.0 if not defined
    if (upfronts.get() == 0) {
        upfronts.reset(new CDoubleArray(n, 0.0));
    }

    // checks all array have the same length n
    if (highStrikes->size() != n || spreads->size() != n || upfronts->size() != n) {
        throw ModelException(method,
            "low strikes, high strikes, spreads and upfronts arrays do not have the same length"); 
    }
    
    // creates new CDOQuotes object with same conventions
    //CDOQuotesSP newQuotes(DYNAMIC_CAST(CDOQuotes, MarketObject::clone()));
    
    // removes market name
    name = "";

    // defines a map from (low strike, high strike) to a maturity ordered set of (maturity, spread, upfront) points
    // so each entry of this map corresponds to a kind of "pseudo tranche quotes" that will be
    // converted into a "real" tranche quotes in a second step below
    typedef map<pair<double, double>, set<TimePoint2D, TimePoint2D::CompareDateC1C2> > StrikesToQuotesMap;
    StrikesToQuotesMap strikesToQuotesMap;
    
    // populates strikesToQuotesMap
    for (i = 0; i < n; ++i) {
		strikesToQuotesMap[pair<double, double>((*lowStrikes)[i], (*highStrikes)[i])].insert(
            TimePoint2D((*expiries)[i]->toDate(valueDate), (*spreads)[i], (*upfronts)[i]));
	}

    // transforms strikesToQuotesMap into CDOTrancheQuotesArray ("pseudo" -> "real" tranche quotes)
    
    int nbTrancheQuotes = strikesToQuotesMap.size();
    CDOTrancheQuotesArraySP newTrancheQuotes(new CDOTrancheQuotesArray(nbTrancheQuotes));

    double tqLowStrike, tqHighStrike; // "tq" stands for tranche quote
    ExpiryArraySP tqExpiries;
    CDoubleArraySP tqSpreads, tqUpfronts;
    int nbTQExpiries;
    
    // loops through "pseudo" tranche quotes
    StrikesToQuotesMap::iterator mapIter = strikesToQuotesMap.begin();
    for (i = 0; mapIter != strikesToQuotesMap.end(); mapIter++){
        tqLowStrike = (*mapIter).first.first;
        tqHighStrike = (*mapIter).first.second;
        nbTQExpiries = (*mapIter).second.size();
        tqExpiries.reset(new ExpiryArray(nbTQExpiries));
        tqSpreads.reset(new CDoubleArray(nbTQExpiries));
        tqUpfronts.reset(new CDoubleArray(nbTQExpiries));
        
        // loops through this "pseudo" tranche quote to populate expiries, spreads and upfronts
        set<TimePoint2D, TimePoint2D::CompareDateC1C2>::iterator setIter = (*mapIter).second.begin();
        for (j = 0; setIter != (*mapIter).second.end(); setIter++) {
			(*tqExpiries)[j].reset(new BenchmarkDate(setIter->getDate()));
            (*tqSpreads)[j] = setIter->getCoord1(); 
            (*tqUpfronts)[j] = setIter->getCoord2(); 

            j++;
		}
        
        // builds "real" tranche quote
        (*newTrancheQuotes)[i] = CDOTrancheQuotes::create(
            valueDate, tqLowStrike, tqHighStrike, tqExpiries, tqSpreads, tqUpfronts);

        i++;
	}
    
    // puts new tranche quotes into newQuotes
    trancheQuotes.reset(new CDOTrancheQuotesWrapperArray(newTrancheQuotes->size()));
    for (i = 0; i < newTrancheQuotes->size(); ++i) {
        // creates wrapper with empty name
        (*trancheQuotes)[i] = CDOTrancheQuotesWrapperSP(new CDOTrancheQuotesWrapper(""));
        
        // populates wrapper
        (*(*trancheQuotes)[i]).setObject((*newTrancheQuotes)[i]);
    }

    // builds strikesToTrancheQuotesIdxMap
    initStrikesToTrancheQuotesIdxMap();
}

/** Main method : builds the CDOquotes object */
CDOQuotesConstSP CDOQuotes::buildCDOQuotes() const
{
    return CDOQuotesConstSP(this);
}



/** Default constructor */
IObject* CDOQuotes::defaultConstructor() {
    return new CDOQuotes();
}

/** TYPE for CDOQuotes */
CClassConstSP const CDOQuotes::TYPE = CClass::registerClassLoadMethod(
    "CDOQuotes", typeid(CDOQuotes), CDOQuotes::load);

DRLIB_END_NAMESPACE
