//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MarketDataFetcher.cpp
//
//   Description : Helper class for models to get data out of market cache
//
//   Author      : Andrew J Swain
//
//   Date        : 1 February 2002
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_MARKETDATAFETCHER_CPP
#include "edginc/MarketData.hpp"
#include "edginc/MarketDataFetcher.hpp"
#include "edginc/Model.hpp"


DRLIB_BEGIN_NAMESPACE

/** Essentially a wrapper around MarketData::GetData(model, name, type).
    Retrieves the market object with
    the specififed name and type (and retrieves its market data). The
    domesticYCName flags what the current 'domestic' currency is (see
    getDomesticYCName() for more information). */
MarketObjectSP MarketDataFetcher::getMarketObject(
    const IModel*     model,
    const MarketData* market,
    const string&     name,
    CClassConstSP     type,
    const string&     domesticYCName) const{
    MarketObjectSP marketObj = model->GetMarket(market, name, type);
    if(marketObj.get()){
        // then get the market data that that object needs
        model->getComponentMarketData(market, marketObj, domesticYCName);
    }
    return marketObj;
}

/** default fetch method. This avoids getting stochastic yield curve unless
 *  configured not to do so. */
MarketObjectSP MarketDataFetcher::fetch(
    const MarketData*    market,
    const string&        name,
    const CClassConstSP& type,
    const IModel*        model) const {
    try {
        // just call standard method to determine type
        CClassConstSP typeToGet = getTypeToRetrieve(type);
        // but need to do this until StochasticYieldCurve type is dropped
        typeToGet = stochasticYCFix(market, name, typeToGet);
        if (!typeToGet){
            return MarketObjectSP();
        }
        return market->GetObjectData(name, typeToGet, model);
    } catch (exception& e) {
        string m("Failed to get market data " + name + " of type " + 
                 type->getName() + (model? " for model of type " +
                                    model->getClass()->getName() : ""));
        throw ModelException(e, "MarketDataFetcher::fetch", m);
    }
}

/** Retrieve a correlation from the cache with the supplied name and
    type specified by corrType. type1 and type2 are the types of the
    two objects with respect to which this is a correlation. The
    return correlation will be correcly configured for the
    appropriate sensitivity (eg PHI, FX_PHI) etc. Note that a model
    may decide that it does not use this type of correlation (eg IR v EQ) 
    and return a 'dummy' version if say the object does not live in the
    cache. In a similar manner, the returned object may be sensitive to
    no sensitivities if the model is not going to use the correlation */
CorrelationBaseSP MarketDataFetcher::getCorrelation(const IModel*     model,
                                                    const string&     corrName,
                                                    CClassConstSP     type1,
                                                    CClassConstSP     type2,
                                                    CClassConstSP     corrType,
                                                    const MarketData* market) const {
    MarketObjectSP corr(market->GetData(model, corrName, corrType));
    if (corr.get()) {
        if (!CorrelationBase::TYPE->isInstance(corr)){
            throw ModelException("MarketDataFetcher::getCorrelation", "Correlation of type "+
                                corrType->getName()+" does not implement "+
                                "CorrelationBase");
        }
        CorrelationBaseSP theCorr(CorrelationBaseSP::dynamicCast(corr));
        theCorr->configureForSensitivities(type1, type2);
        model->getComponentMarketData(market, corr);
        return theCorr;
    }
    else {
        return CorrelationBaseSP(0);
    }
}

/** Invoked for each piece of market data (whether already inline in
    instrument or pulled from cache). The default implementation just
    returns mo. Derived classes can use this to replace market data
    objects with other instances. It compliments GetMarket in that this
    method works for instruments which have the market data inside them
    already */
MarketObjectSP MarketDataFetcher::modifyMarketData(
    const IModel*         model,
    const MarketData*     market,
    const CClassConstSP&  clazz,     // type originally requested
    const MarketObjectSP& mo) const  /* what GetMarket returned or
                                      * what was "inline" already */
{
    return mo;
}

//// throws an exception if all types have been loaded
void checkThatTypeRegistrationIsInProgress(){
    if (!CClass::classLoadingInProgress()){
        throw ModelException("MarketDataFetcher::setDefaultRetrievalMode",
                             "This method may only be called during Classes"
                             " load method");
    }
}

/** Class to capture the configuration of MarketDataFetchers. It captures
    what types to get, when to get them etc. */
class MarketDataFetcher::RetrievalConfig{
public: // visible to this source file
    CClassConstSP requestedType;  // eg IRVol
    CClassConstSP containingType; // eg YieldCurve, can be null
    bool          retrieve;       // false=>don't get
    CClassConstSP typeToGet;      /* if retrieve is true then specialised type,
                                     use NULL to use actual type requested in
                                     the fetch method */

    // for use with hash_set
    struct HashAndEqual{
        // hash function for key
        size_t operator()(const RetrievalConfig& config) const {
            return (size_t)config.requestedType;
        }
        // equals function for key
        bool operator()(const RetrievalConfig& config1,
                        const RetrievalConfig& config2) const {
            return (config1.requestedType == config2.requestedType);
        }
    };
    /** Given an object of type 'type' which implements/is derived from both
        parentType1 and parentType2 then determine which parent is the
        'closer' parent. This basically means determining if parentType1 is
        derived from parentType2 (or vice versa) and then returning the
        'child' type. Note that if parentType1 and parentType2 are not
        related (ie neither derived from each other) then an exception is
        thrown */
    static CClassConstSP resolveTypeAmbiguity(CClassConstSP type,
                                              CClassConstSP parentType1,
                                              CClassConstSP parentType2){
        // possible ambiguity - make sure it's ok
        bool useParent1 = parentType2->isAssignableFrom(parentType1);
        bool useParent2 = parentType1->isAssignableFrom(parentType2);
        if (useParent1 ^ useParent2){
            return useParent1? parentType1: parentType2;
        }
        string m("Ambiguity when processing Market Data retrieval configuration"
                 "for type "+type->getName()+
                 ". Configuration information specified for types "+
                 parentType1->getName()+" and "+parentType2->getName());
        throw ModelException("MarketDataFetcher::RetrievalConfig::"
                             "resolveTypeAmbiguity", m);
    }

    /** for each requested type may need several configuration entries
        for being inside different objects. So use a vector<vector<>>.
        Note no real advantage of using a set */
    typedef vector<vector<RetrievalConfig> > Table;

    RetrievalConfig(CClassConstSP requestedType,
                    CClassConstSP containingType,
                    bool          retrieve,
                    CClassConstSP typeToGet,
                    bool          validateTypes):
        requestedType(requestedType), containingType(containingType),
        retrieve(retrieve), typeToGet(typeToGet){
        static const string method("MarketDataFetcher::RetrievalConfig");
        if (!retrieve && typeToGet){
            throw ModelException(method, "Type to get cannot be specifed if "
                                 "retrieve flag is set to false");
        }
        if (validateTypes && typeToGet && 
            !requestedType->isAssignableFrom(typeToGet)){
            string m("Can't specialise requests for type "+
                     requestedType->getName()+" with type "+
                     typeToGet->getName());
            throw ModelException(method, m);
        }
    }

    /** Finds most relevant entry in a set of RetrievalConfigs for the 
        specified type. Fails if there is ambiguity (eg class A implements
        interfaces I1 and I2, but I1 and I2 are not connected and there are
        different rules for I1 and I2). Returns null if there is no entry */
    //static pair<HashSet::iterator, HashSet::iterator> findEntry(
    static const vector<RetrievalConfig>* findEntry(
        const Table&               retrievalConfigTable,
        CClassConstSP              requestedType)
    {
        static const string method("MarketDataFetcher::RetrievalConfig::"
                                   "findEntry");
        const vector<RetrievalConfig>* currentEntry = NULL;
        for (unsigned int i = 0; i < retrievalConfigTable.size(); i++){
            if (retrievalConfigTable[i].size() > 0){
                CClassConstSP entryType = 
                    retrievalConfigTable[i][0].requestedType;
                if (entryType->isAssignableFrom(requestedType)){
                    // found a suitable Config entry
                    if (currentEntry){
                        // possible ambiguity - make sure it's ok
                        CClassConstSP currentEntryType = 
                            (*currentEntry)[0].requestedType;
                        if (resolveTypeAmbiguity(requestedType,
                                                 currentEntryType,
                                                 entryType) == entryType){
                            currentEntry = &retrievalConfigTable[i];
                        }
                    } else {
                        currentEntry = &retrievalConfigTable[i];
                    }
                }
            }
        }
        return currentEntry;
    }

    //// add the supplied retrievalConfig into the retrievalConfigTable.
    //// Existing entries are overrwritten as appropriate
    static void add(Table&                  retrievalConfigTable,
                    const RetrievalConfig&  retrievalConfig){
        CClassConstSP entryType = retrievalConfig.requestedType;
        // search for entry Type
        for (unsigned int i = 0; i < retrievalConfigTable.size(); i++){
            vector<RetrievalConfig>& currentEntry = retrievalConfigTable[i];
            if (currentEntry.size() > 0){
                // note do an exact match on type
                if (currentEntry[0].requestedType == entryType){
                    // now see if containingType matches
                    for (unsigned int j = 0; j < currentEntry.size(); j++){
                        if (retrievalConfig.containingType == 
                            currentEntry[j].containingType){
                            // then just overwrite retrieve and typeToGet
                            currentEntry[j].retrieve = retrievalConfig.retrieve;
                            currentEntry[j].typeToGet = 
                                retrievalConfig.typeToGet;
                            return;
                        }
                    }
                    currentEntry.push_back(retrievalConfig);
                    return;
                }
            }
        }
        retrievalConfigTable.
            push_back(vector<RetrievalConfig>(1, retrievalConfig));
    }

    /** Find the entry for the specified type (NB matching exactly on type).
        This is the partner to the add method above. Returns NULL if
        no entry found */
    static const RetrievalConfig* get(const Table&  retrievalConfigTable,
                                      CClassConstSP entryType,
                                      CClassConstSP containingType){
        // search for entry Type
        for (unsigned int i = 0; i < retrievalConfigTable.size(); i++){
            const vector<RetrievalConfig>& currentEntry =
                retrievalConfigTable[i];
            if (currentEntry.size() > 0){
                // note do an exact match on type
                if (currentEntry[0].requestedType == entryType){
                    // now see if containingType matches
                    for (unsigned int j = 0; j < currentEntry.size(); j++){
                        if (containingType == currentEntry[j].containingType){
                            return &currentEntry[j];
                        }
                    }
                    return NULL;
                }
            }
        }
        return NULL;
    }

    /** Given a RetrievalConfig and the original requestedType what do
        we get. NULL implies don't get */
    CClassConstSP getTypeToRetrieve(CClassConstSP requestedType) const{
        if (!retrieve){
            return NULL;
        }
        return typeToGet? typeToGet: requestedType;
    }

    /** Given the requested type, and a corresponding RetrievalConfig,
        works out what type to actually get. A null return value means
        not to retrieve the item at all */
    static CClassConstSP getTypeToRetrieve(
        CClassConstSP                  requestedType,
        const vector<RetrievalConfig>* entry,
        const vector<CClassConstSP>&   objectTypeStack)
    {
        if (!entry || entry->empty()){
            // simple case - just get what was asked for
            return requestedType;
        }
        if (entry->size() == 1 && !(*entry)[0].containingType){
            // simple too - entry independent of objectTypeStack
            CClassConstSP typeToGet = entry->front().typeToGet;
            if (!typeToGet || requestedType->isAssignableFrom(typeToGet)){
                // should only consider case where requested type is derived
                // from the override type (ie typeToGet) if not null
                return entry->front().getTypeToRetrieve(requestedType);
            } else {
                // no rules apply, so just get what was asked for
                return requestedType;
            }
        }
        if (objectTypeStack.empty()){
            // if we're not inside any component then we should just use
            // the rule, if specified, where the containingType is NULL
            for (unsigned int i = 0; i < entry->size(); i++){
                if (!(*entry)[i].containingType){ // found default entry
                    return (*entry)[i].getTypeToRetrieve(requestedType);
                }
            }
            return requestedType; // no rule
        }

        // now it gets more complicated. We have to search through the
        // objectTypeStack to see if any objects match up with the
        // containingType inside the vector<RetrievalConfig>. Loop
        // backwards since the higher index corresponds to the top of
        // stack
        const RetrievalConfig* defaultConfig = 0; /* will hold RetrievalConfig
                                                     with null containingType */
        const RetrievalConfig* currentEntry = 0; /* will hold RetrievalConfig 
                                                    to use */
        /* Note: order of loop through objectTypeStack is from outer most object
           first as this seems more useful in practice (eg don't get IR vols in
           yield curves unless the yield curve is in a QuantoCDSParSpreads 
           object). If more flexibility is needed, consider use of 'weak' or
           'strong' designation for an entry. If weak then keep on looping
           through objectTypeStack rather than returning early */
        for (unsigned int i = 0; i < objectTypeStack.size(); i++){
            CClassConstSP objectType = objectTypeStack[i];
            for (unsigned int j = 0; j < entry->size(); j++){
                // should only consider case where requested type is derived
                // from the override type (ie typeToGet) if not null
                CClassConstSP typeToGet = (*entry)[j].typeToGet;
                if (!typeToGet || requestedType->isAssignableFrom(typeToGet)){
                    CClassConstSP configContainingType = 
                        (*entry)[j].containingType;
                    if (!configContainingType){
                        defaultConfig = &(*entry)[j];
                    } else if (configContainingType->
                               isAssignableFrom(objectType)){
                        if (!currentEntry){
                            currentEntry = &(*entry)[j];
                        } else {
                            // have possible ambiguity
                            if (resolveTypeAmbiguity(
                                    objectType,
                                    currentEntry->containingType,
                                    configContainingType) == 
                                configContainingType){
                                currentEntry = &(*entry)[j];
                            }
                        }
                    }
                }
            }
            if (currentEntry){
                // found a match - use it
                return currentEntry->getTypeToRetrieve(requestedType);
            }
        }
        if (defaultConfig){
            // no relevant containing type - use default value for this type
            return defaultConfig->getTypeToRetrieve(requestedType);
        }
        // no information for this type. Return original type requested
        return requestedType;
    }
};

/** Implementation details for MarketDataFetcher - allows us to change fields/
    methods etc without changing the header file */
class MarketDataFetcher::Imp{
public:
    // the static Table is what is used to initialise each MarketDataFetcher
    static RetrievalConfig::Table defaultRetrievalConfig;
    RetrievalConfig::Table        retrievalConfig;
    // as we recursively retrieve market data, this the stack of objects
    // that we are 'inside'
    vector<CClassConstSP>         objectTypeStack;
    /* what the current 'domestic' yield curve is, if available, eg
       the instrument discount curve */
    string                        domesticYCName;
    // work out what to do with this
    DateTimeSP                    corrSwapExpiry; // default: false

    //// constructor copies defaultRetrievalConfig into its own copy
    Imp():retrievalConfig(defaultRetrievalConfig){}
};

// definition of static field
MarketDataFetcher::RetrievalConfig::Table 
MarketDataFetcher::Imp::defaultRetrievalConfig;

//// this vanishes once get rid of StochasticYieldCurves
CClassConstSP MarketDataFetcher::stochasticYCFix(
    const MarketData*    market,
    const string&        name,
    CClassConstSP        type) const {

    if (!type){
        return type; // cope with NULL type
    }
    CClassConstSP modType = type;
    static CClassConstSP iYCType = CClass::forName("IYieldCurve");
    static CClassConstSP iStochYCType =
        CClass::forName("IStochasticYieldCurve");
    static CClassConstSP iDetermYCType =
        CClass::forName("IDeterministicYieldCurve");
    static CClassConstSP irVolBaseType = CClass::forName("IRVolBase");

    if (iYCType->isAssignableFrom(type)){
        // if we were in a yield curve, do we get the irVol? A bit hacky :-)
        my->objectTypeStack.push_back(iYCType);
        bool getStochasticYCs = getTypeToRetrieve(irVolBaseType) != NULL;
        my->objectTypeStack.pop_back();
        
        if (!getStochasticYCs){
            // been asked for a yield curve that isn't a stochastic one
            modType = iDetermYCType; // get a determinstic one
        } else if (getStochasticYCs && 
                   !iDetermYCType->isAssignableFrom(type)){
            // been asked for a yield curve that isn't a determinstic one.
            // Is there a IStochasticYieldCurve in the market data cache?
            if (market->hasData(name, iStochYCType)){
                // yes so use it
                modType = iStochYCType; // get a stochastic one
            } else {
                // just get an ordinary yield curve which we expect to
                // have the irVol in it
            }
        }
    }
    return modType;
}

/** Select whether market data of the specified type should be, by
    default, retrieved or not from the MarketData cache and, if it is,
    what specialised type should be retrieved. Use NULL for typeToGet
    to indicate that the objects of requestedType should have the type
    requested unchanged. This is a global setting and this method can
    only be called at startup (ie during load methods).  It
    [potentially] affects all market data fetchers so care needs to be
    exercised (otherwise you may break every test case ...) */
void MarketDataFetcher::setDefaultRetrievalMode(CClassConstSP requestedType,
                                                bool          retrieve,
                                                CClassConstSP typeToGet){
    checkThatTypeRegistrationIsInProgress();
    // don't validate type info, since not all types loaded at this point
    RetrievalConfig::add(Imp::defaultRetrievalConfig,
                         RetrievalConfig(requestedType, 0, retrieve,
                                         typeToGet, false));
}

/** This is the same as the above method but the rule for retrieving or
    not retrieving the object is only used if, when retrieving data, we are
    inside an object of the specified containing type. For example, this is
    used to avoid getting IRVols inside yield curves. This is a global setting
    and this method can only be called at startup (ie during load methods).
    It [potentially] affects all market data fetchers so care needs to be
    exercised (otherwise you may break every test case ...)  */
void MarketDataFetcher::setDefaultRetrievalMode(CClassConstSP requestedType,
                                                CClassConstSP containingType,
                                                bool          retrieve,
                                                CClassConstSP typeToGet){
    checkThatTypeRegistrationIsInProgress();
    // don't validate type info, since not all types loaded at this point
    RetrievalConfig::add(Imp::defaultRetrievalConfig,
                         RetrievalConfig(requestedType, containingType,
                                         retrieve, typeToGet, false));
}


/** Select whether market data of the specified type should be
    retrieved or not from the MarketData cache by this
    MarketDataFetcher. This is the same as setDefaultRetrievalMode but is
    specific to this instance of the MarketDataFetcher  */
void MarketDataFetcher::setRetrievalMode(
    CClassConstSP requestedType, 
    bool          retrieve,
    CClassConstSP typeToGet) const /* hacky */
{
    RetrievalConfig::add(my->retrievalConfig, /* with the amazing powers of C++
                                                 this is not const! */
                         RetrievalConfig(requestedType, 0, 
                                         retrieve, typeToGet, true));
}

/** Query whether market data of the specified type is retrieved or
    not from the MarketData cache by this MarketDataFetcher ignoring
    any settings dependent on what object we are in when retrieving
    market data. This is the partner to setRetrievalMode method which
    takes 3 arguments. */
bool MarketDataFetcher::getRetrievalMode(CClassConstSP  requestedType) const {
    const RetrievalConfig* config = RetrievalConfig::get(my->retrievalConfig, 
                                                         requestedType, NULL);
    return !config? true: config->retrieve;
}

/** Query whether market data of the specified type is retrieved or
    not from the MarketData cache by this MarketDataFetcher when
    inside an object of type containingType. This is the partner to
    setRetrievalMode method which takes 4 arguments. */
bool MarketDataFetcher::getRetrievalMode(CClassConstSP  requestedType,
                                         CClassConstSP  containingType) const{
    const RetrievalConfig* config = RetrievalConfig::get(my->retrievalConfig, 
                                                         requestedType, 
                                                         containingType);
    return !config? true: config->retrieve;
}

/** This is the same as the above method but the rule for retrieving or
    not retrieving the object is only used if, when retrieving data, we are
    inside an object of the specified containing type. For example, this is
    used to avoid getting IRVols inside yield curves  */
void MarketDataFetcher::setRetrievalMode(
    CClassConstSP requestedType,
    CClassConstSP containingType,
    bool          retrieve,
    CClassConstSP typeToGet) const /* hacky */
{
    RetrievalConfig::add(my->retrievalConfig,
                         RetrievalConfig(requestedType, containingType, 
                                         retrieve, typeToGet, true));
}

/** Given the requested type, works out what type to actually get. A null
    return value means not to retrieve the item at all */
CClassConstSP MarketDataFetcher::getTypeToRetrieve(
    CClassConstSP requestedType) const{
    const vector<RetrievalConfig>* entry = 
        RetrievalConfig::findEntry(my->retrievalConfig, requestedType);
    return RetrievalConfig::getTypeToRetrieve(requestedType, entry,
                                              my->objectTypeStack);
}

MarketDataFetcher::~MarketDataFetcher(){}

MarketDataFetcher::MarketDataFetcher(): my(new Imp()){}

/** Resolve an array of multiply selected market data
    into a single chosen market datum */
MarketObjectConstSP MarketDataFetcher::resolveMultiData(
    const MarketObjectArray& marketObjs,
    CClassConstSP            type) const
{
    static const string routine = "MarketDataFetcher::resolveMultiData";

    string typeName = type->getName();
    string rootName;
    MarketObjectArray exactMarketObjs, derivedMarketObjs;
    MarketObjectConstSP marketObj;

    if (!marketObjs.empty()) rootName = marketObjs[0]->getRootName();

    for (MarketObjectArray::const_iterator i = marketObjs.begin();
        i != marketObjs.end(); ++i)
    {
        if ((*i)->getClass()->getName() == typeName)
            exactMarketObjs.push_back(*i);
        else
            derivedMarketObjs.push_back(*i);
    }

    if (exactMarketObjs.size() > 1) {
        // trouble - there are at least 2 exact matching types
        string s("Object with name " + rootName + " and type " +
                    typeName + " is not unique in cache. "
                    "At least two objects with type " +
                    typeName + " exist" );
        throw ModelException(routine, s);
    }
    else if (exactMarketObjs.size() == 1) {
        marketObj = exactMarketObjs[0];
    }
    else if (derivedMarketObjs.size() > 1) {
        // trouble - there are at least 2 matching types but none of them exact
        string s("Object with name " + rootName + " and type " +
                typeName + " is not unique in cache. "
                "At least two objects with types " +
                derivedMarketObjs[0]->getClass()->getName() + " and " +
                derivedMarketObjs[1]->getClass()->getName() + " exist" );
        throw ModelException(routine, s);
    }
    else if (derivedMarketObjs.size() == 1) {
        marketObj = derivedMarketObjs[0];
    }

    return marketObj;
}

/** By default, the same method on Model delegates to this method.
    Briefly, it allows the MarketDataFetcher context information 
    regarding what object it is currently in */
void MarketDataFetcher::getComponentMarketData(const IModel*        model,
                                               const MarketData*    market,
                                               MarketObjectSP       mo) const
{
    if (!mo.get())
    {
        throw ModelException("MarketDataFetcher::getComponentMarketData",
                             "Could not get data for MarketObject!");
    }
    my->objectTypeStack.push_back(mo->getClass());
    mo->getMarket(model, market);
    my->objectTypeStack.pop_back();
}


// need to work out what to do with this
DateTimeSP MarketDataFetcher::setCorrSwapExpiry(DateTimeSP corrSwapExpiryInput) const {
    DateTimeSP oldValue = my->corrSwapExpiry;
    my->corrSwapExpiry = corrSwapExpiryInput;
    bool retrieve = corrSwapExpiryInput.get() != NULL;
    setRetrievalMode(CClass::forName("CorrSwapSamplingAdj"), retrieve, NULL);
    setRetrievalMode(CClass::forName("CorrSwapBasisAdj"), retrieve, NULL);
    setRetrievalMode(CClass::forName("CorrelationCategory"), retrieve, NULL);
    return oldValue;
}
DateTimeSP MarketDataFetcher::getCorrSwapExpiry() const {
    // need to work out what to do here - value is used in Correlation.cpp
    return my->corrSwapExpiry;
}

/** Same as above but specifies a 'domestic' yield curve for the 
    current component (see getDomesticYCName()) */
/* virtual */ void MarketDataFetcher::getComponentMarketData(
    const IModel*     model,
    const MarketData* market,
    MarketObjectSP    mo,
    const string&     domesticYCName) const{
    string previousYCName = my->domesticYCName; // save
    my->domesticYCName = domesticYCName; // set
    getComponentMarketData(model, market, mo);
    my->domesticYCName = previousYCName; // restore
}

/** Return the name of the 'domestic' yield curve, if available, or
 * an empty string if not available. */
const string& MarketDataFetcher::getDomesticYCName() const{
    return my->domesticYCName;
}

/** Sets the 'domestic' yield curve name */
void MarketDataFetcher::setDomesticYCName(string domYCName) const {
    my->domesticYCName = domYCName;
}

DRLIB_END_NAMESPACE
