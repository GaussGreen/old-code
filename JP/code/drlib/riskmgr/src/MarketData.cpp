//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MarketData.cpp
//
//   Description : Class controlling access to market data
//
//   Author      : Mark A Robson
//
//   Date        : 16 Jan 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_MARKETDATA_CPP
#include "edginc/MarketData.hpp"
//#include "edginc/Model.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Writer.hpp"
#include "edginc/Hashtable.hpp"
#include "edginc/ObjectIteration.hpp"
#include "edginc/MarketDataFetcher.hpp"
#include "edginc/NonPricingModel.hpp"
#include "edginc/IMarketObjectQualifier.hpp"
#include "edginc/ClientRunnable.hpp"
#include <map>
#include <fstream>
#include ext_hash_map

DRLIB_BEGIN_NAMESPACE

///////// MarketData class //////////

class MarketData::Imp{
public:
    ~Imp(){}

    // simple constructors
    Imp(){}

    Imp(const DateTime& today): RefDate(today){}
    
    /** Class for identifying something by a pair of string - order of which
        may or may not be irrelevant */
    class StringPair{

    public:
        /** constructor */
        StringPair(const string& name1, const string& name2):
            id1(name1), id2(name2){}
        /** Symmetric hash functions */
        struct SymHashFuncs{
            /** equality operator */
            bool operator()(const StringPair& name1, 
                            const StringPair& name2) const{
                return ((name1.id1 == name2.id1 && name1.id2 == name2.id2) ||
                        (name1.id2 == name2.id1 && name1.id1 == name2.id2));
            }
            
            /** hash operator */
            size_t operator()(const StringPair& name) const{
                // hash needs to be symmetric under id1 and id2
                return hash_string(name.id1.c_str()) ^ 
                    hash_string(name.id2.c_str());
            }
        };
        /** Asymmetric hash functions */
        struct ASymHashFuncs{
            /** equality operator */
            bool operator()(const StringPair& name1, 
                            const StringPair& name2) const{
                return (name1.id1 == name2.id1 && name1.id2 == name2.id2);
            }
            /** hash operator */
            size_t operator()(const StringPair& name) const{
                return hash_string(name.id1.c_str()) ^
                    hash_string(name.id2.c_str());
            }
        };
        string id1;
        string id2;
    };

    typedef multimap<const string, MarketObjectSP> CDataMap;

    typedef hash_map<StringPair, string, 
        StringPair::SymHashFuncs, StringPair::SymHashFuncs> SymStringPairHash;
    // We may want to parameterise the type of data we're storing eg
    // pass in "Correlation" etc rather than have a separate function for
    // each one

    typedef hash_map<StringPair, string, StringPair::ASymHashFuncs, 
        StringPair::ASymHashFuncs> ASymStringPairHash;

    typedef hash_map<string, string, Hashtable::StringHash> StringTable;

    DateTime RefDate;
    /* Data identified by key together with its type ie 'duplicate' entries in
       the multimap correspond to objects with the same name but different
       type */
    CDataMap DataCache;
    
    /** Allows us to identify correlations by supplying the names of the
        two objects which are correlated */
    SymStringPairHash correlationsHash;
    SymStringPairHash correlationTermHash;
	SymStringPairHash corrSwapSamplingAdjHash;
    /** Allows us to resolve the name of the fx asset which relates these
        two yield curves */
    ASymStringPairHash fxHash;
    StringTable isoCodesByYC; // iso code for each YC

    /**
     * App-specific caches: see MarketData::plugin()
     *
     * Left NULL on construction, and NULLed by MarketData::AddData to force
     * re-cacheing whenever our contents change.  MarketData::clone() sets the
     * cloned MarketData to share plugins initially (which is sound because
     * plugins are supposed to be computed purely from the MarketData contents
     * which are immutable, and because of the above-noted clear-on-change).
     *
     * See MarketDataPluginCacheTest below.
     */

    refCountPtr<hash_map<string, IObjectConstSP, Hashtable::StringHash> > plugins;

    /* Data for ObservableHistorys - identified by key (name of MarketObservable)
       together with its source and type NB 'duplicate' entries in
       the multimap correspond to histories for the same observable but with 
       different source and we don't allow the same observable to have 
       histories of different types*/
    class ObsHistoryDetails {
    public: 
        ObsHistoryDetails(const string& source, const string& cacheName,
            const CClassConstSP& histType) : source(source), 
                      cacheName(cacheName), histType(histType) {}

        string          source;
        string          cacheName;
        CClassConstSP   histType;
    };
    typedef refCountPtr<ObsHistoryDetails> ObsHistoryDetailsSP;

    typedef multimap<const string, ObsHistoryDetailsSP> ObsHistMap;
    ObsHistMap obsHistory;
};

// static fields //
const string MarketData::REF_DATE = "RefDate";
const string MarketData::DATA_CACHE = "DataCache";

MarketData::~MarketData(){}

MarketData::MarketData(): CObject(TYPE), my(new Imp()){}

MarketData::MarketData(const DateTime &date): 
    CObject(TYPE), my(new Imp(date)){}

/** creates copy of MarketData but does not clone components */
IObject* MarketData::clone() const{
    MarketData* clone = 0;
    try{
        clone = new MarketData(my->RefDate);
        // then iterate through copying over entries. Note (1) no deep
        // copy. (2) Important to go through AddData method to ensure
        // other fields are set up
        for (Imp::CDataMap::const_iterator myIterator = my->DataCache.begin();
             !(myIterator == my->DataCache.end());
             ++myIterator){
            clone->AddData(myIterator->second);
        }

        // Ensure that the two MarketData's share a cache, so that each can
        // use caches created for the other, until one or the other's
        // contents change --- see comments on Imp::plugins above

        if (!my->plugins) {
            // need to ensure that we really are sharing plugins even if
            // none so far
            my->plugins.reset(new hash_map<string, IObjectConstSP,
                                           Hashtable::StringHash>());
        }

        clone->my->plugins = my->plugins;
    } catch (exception& e){
        throw ModelException(e, "MarketData::clone",
                             "Failed to clone MarketData");
    }
    return clone;
}

/** write object out to writer */
void MarketData::write(const string& tag, Writer* writer) const {
    try {
        IObjectConstSP obj(writer->objectStart(tag, "", this, true));
        if (obj.get()){
            // if we have implemented IPrivateObject then we need to act
            // on obj here which will be of type IPublicObject
            my->RefDate.write(REF_DATE, writer);
            writer->objectStart(DATA_CACHE, "", this, true);
            for (Imp::CDataMap::const_iterator iter = my->DataCache.begin();
                 !(iter == my->DataCache.end()); ++iter){
                iter->second->write("entry", writer);
            }
            writer->objectEnd(DATA_CACHE, this);
        }
        writer->objectEnd(tag, this);
    }
    catch (exception& e){
        throw ModelException(&e, "MarketData::write");
    }
}

/** populate an empty object from reader */
void MarketData::import(Reader::Node* elem, Reader* reader){
    const static string routine = "MarketData::import";
    try{
        Reader::NodeListSP nl(elem->children());
        for (unsigned int i = 0; i < nl->size(); i++) {
            Reader::Node* child = (*nl)[i].get();
            string fieldName = child->name();
            if (fieldName == REF_DATE){
                IObjectSP obj(reader->read(child));
                try{
                    my->RefDate = *(DateTimeSP::dynamicCast(obj));
                } catch (exception& e){
                    throw ModelException(e, "Failed to read date from"
                                         " object of type "+
                                         obj->getClass()->getName());
                }
            } else if (fieldName == DATA_CACHE){
                Reader::NodeListSP mktData(child->children());
                for (unsigned int i = 0; i < mktData->size(); i++) {
                    Reader::Node* mktChild = (*mktData)[i].get();
                    // read in entry
                    IObjectSP object(reader->read(mktChild));
                    try{
                        AddData(MarketObjectSP::dynamicCast(object));
                    } catch (exception& e){
                        string m("Failed to add object of type "+
                                 object->getClass()->getName());
                        throw ModelException(e, m);
                    }
                }
            }
        }
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}

/** write object out in 'output' format - ie suitable for comparing
    regression files with */
void MarketData::outputWrite(const string& linePrefix,
                             const string& prefix, ostream& stream) const{
    my->RefDate.outputWrite(linePrefix, 
                            prefix.empty()? REF_DATE: prefix + REF_DATE,
                            stream);
    for (Imp::CDataMap::const_iterator iter = my->DataCache.begin();
         !(iter == my->DataCache.end()); ++iter)
    {
        iter->second->outputWrite(linePrefix,
                                  prefix.empty()? 
                                  iter->first: prefix + "_" + iter->first,
                                  stream);
    }
}

/** Stores (a reference) to the supplied data in this market data
    cache. The getRootName method() on data is used for the object's id */
void MarketData::AddData(const MarketObjectSP& data)
{
    try {
        data->initialise(this);
        string name = data->getRootName();
        IMarketObjectQualifierSP moq1 = data->getMarketObjectQualifier();
        typedef Imp::CDataMap::iterator I;
        pair<I,I> limits = my->DataCache.equal_range(name);
        bool found = false;
        for (I iter = limits.first; !(iter == limits.second) && !found; ++iter){
            if (iter->second->getClass() == data->getClass()){
                IMarketObjectQualifierSP moq2 = iter->second->getMarketObjectQualifier();
                if (!moq1.get() && !moq2.get() ||
                    moq1.get() && moq2.get() && moq1->equals(moq2.get())) {
                    // entry already exists - overwrite
                    iter->second = data;
                    found = true;
                }
            }
        }
        if (!found){
            // insert new entry
            /* in theory can just do "DataCache.insert(make_pair(name, data));"
               but doesn't compile on MSVC. */
            Imp::CDataMap::value_type entry(name, data);
            my->DataCache.insert(entry);
        }

        my->plugins.reset();

    } catch (exception& e){
        throw ModelException(e, "MarketData::AddData");
    }
}
/** Merge supplied market with this one. Data in market overwrites data in
    this if conflicts exist */
void MarketData::addData(const MarketData& market){
    try{
        typedef Imp::CDataMap::const_iterator I;
        for(I iter(market.my->DataCache.begin()); 
            iter != market.my->DataCache.end();
            ++iter){
            AddData(iter->second);
        }
    } catch (exception& e){
        throw ModelException(e, "MarketData::addData");
    }
}


MarketObjectArraySP MarketData::getEquityFXCorrelations(
    const string& equityName,
    const CClassConstSP&  fxAssetType,
    const CClassConstSP&  correlationType)
{
    MarketObjectArraySP marketObjArray(new MarketObjectArray(0));
    string secondAsset;

    for (Imp::SymStringPairHash::const_iterator myIterator =
             my->correlationsHash.begin();
         !(myIterator == my->correlationsHash.end()); ++myIterator) {
        Imp::StringPair hashKey = myIterator->first;

        if ( hashKey.id1 == equityName) {
            secondAsset = hashKey.id2;
        } else if ( hashKey.id2 == equityName ) {
            secondAsset = hashKey.id1;
        }

        /* that's well inefficient, but MSVC complains about ambiguity
           in the constructor when I use DataCache.equal_range as in
           MarketData::getData - will have to revisit this */
        if ( hashKey.id1 == equityName || hashKey.id2 == equityName ) {
            if ( hasData(secondAsset, fxAssetType)) {
                /* the correlation object is representing an equity-fx
                   correlation */
                MarketObjectSP clonedObj = GetData(myIterator->second,
                                                   correlationType);
                marketObjArray->push_back(clonedObj);
            }
        }
    }
    return marketObjArray;
}

/** Returns a [deep] copy of all objects in the cache with the 
    the supplied name. */
MarketObjectArraySP MarketData::GetAllDataWithName(const string& name)const{
    static const string routine("MarketData::GetAllDataWithName");
    try{
        MarketObjectArraySP marketObjects(new MarketObjectArray(0));
        typedef Imp::CDataMap::const_iterator I;
        // search all entries with the same name to see if we can get a match
        pair<I,I> limits = my->DataCache.equal_range(name);
        for (I iter = limits.first; !(iter == limits.second); ++iter){
            IObjectSP clone(iter->second->clone());
            clone = CObject::convertToPublicRep(clone);
            MarketObjectSP clonedObj(MarketObjectSP::dynamicCast(clone));
            marketObjects->push_back(clonedObj);
        }

        return marketObjects;
    } catch (exception& e){
        throw ModelException(e, routine, "Failed to retrieve objects with "
                             "name "+name);
    }
}

IObjectConstSP MarketData::_plugin(
        const type_info& which,
        IObjectConstSP (*factory)(const MarketData&)) const {
    if (!my->plugins) {
        my->plugins.reset(new hash_map<string, IObjectConstSP,
                                       Hashtable::StringHash>());
    }

    IObjectConstSP& it = (*my->plugins)[which.name()];
    if (!it) it = (*factory)(*this);
    return it;
}

/** Returns a [shallow] copy of all objects in the cache with the 
    the supplied type. */
MarketObjectArraySP MarketData::GetAllDataWithType(const CClassConstSP& type)const{
    static const string routine("MarketData::GetAllDataWithName");
    try{
        MarketObjectArraySP marketObjects(new MarketObjectArray(0));
        typedef Imp::CDataMap::const_iterator I;
        // search all entries  to see if we can get a match
        MarketObjectSP obj;
        for (I iter = my->DataCache.begin(); iter != my->DataCache.end(); ++iter) {
			if (type->isInstance(iter->second.get())){
				obj = MarketObjectSP((iter->second));
				marketObjects->push_back(obj);
			}
        }

        return marketObjects;
    } catch (exception& e) {
        throw ModelException(e, routine, "Failed to retrieve objects with "
                             "type "+ type->getName());
    }
}

/** Returns a reference to the requested MarketObject object from
    this cache.  If there is a unique object with the given root name and
    of the given type (or derived from the given type) then that object is
    returned.  The model has the opportunity to resolve multiple objects to
    a unique object.  Otherwise an exception is thrown. */
MarketObjectConstSP MarketData::getData(const string&         rootName,
                                        const CClassConstSP&  type,
                                        const IModel*         model) const{
    static const string routine("MarketData::GetData"); // for backwards compat
    try{
        typedef Imp::CDataMap::const_iterator I;
        // search all entries with the same name
        // amd try to resolve to a unique match
        MarketObjectArray marketObjs;
        MarketObjectConstSP marketObj;
        string typeName = type->getName();
#define CHECK_CACHE 0
#if CHECK_CACHE
        std::ofstream out("c:/temp/marketData.txt");
        out << "Looking for root name: " << rootName << "/ttype: " << type->getName() << "\n";
        I iii= my->DataCache.begin();
        for ( ; iii != my->DataCache.end(); ++iii ) {
            MarketObjectSP candidate(iii->second);
            out << "entry name: " << iii->first
                << "\t entry type: " << candidate->getClass()->getName() << "\n";
        }
#endif
        pair<I,I> limits = my->DataCache.equal_range(rootName);
#if CHECK_CACHE
        out << "got candidates\n";
#endif
        for (I iter = limits.first; !(iter == limits.second); ++iter){
            MarketObjectSP candidate(iter->second);
#if CHECK_CACHE
            out << "candidate name: " << rootName
                << "\t candidate type: " << candidate->getClass()->getName() << "\n";

#endif
            /* we support some cases of multiple objects with the same name
               and both derived from supplied type (eg FTSE surface and FTSE 
               parameterised vol, when asking for LN vol)

               in the absence of a model (model == NULL) we allow:
                    - one and only exact or derived type object: return it
                    - one exact and many derived type objects: return exact type object
                    thus we don't allow:
                    - multiple exact type objects
                    - multiple derived type objects only

                in the presence of a model:
                    - model resolves, return result
            */
            if (type->isInstance(candidate.get())){
                // found a promising candidate

                // does the object use a proxy (eg VolPreferred)
                CClassConstSP proxyType = candidate->proxyType(this, model);
                if (proxyType){
                    if (proxyType->isAssignableFrom(candidate->getClass())){
                        // avoid infinite recursion
                        string m(rootName+"'s request to use an "
                                 "object of type "+proxyType->getName()+
                                 " is not specific enough");
                        throw ModelException(routine, m);
                    }
                    // replace with proxy
                    candidate = MarketObjectSP::constCast(
                        getData(rootName, proxyType)); // recursion
                }

                marketObjs.push_back(candidate);
            }
        }
#if CHECK_CACHE
        out.close();
#endif
        marketObj = model != NULL ?
            model->getMDF()->resolveMultiData(marketObjs, type) :
            NonPricingModel().getMDF()->resolveMultiData(marketObjs, type);

        if (!marketObj){
            throw ModelException(routine, "Object not found in cache");
        }
        return marketObj;
    } catch (exception& e){
        throw ModelException(e, routine, "Failed to retrieve object with "
                             "name "+rootName+" and type "+type->getName());
    }
}

/** Returns a [deep] copy of the requested MarketObject object from
    this cache. If a unique object with the given name and of the
    given type (or derived from the given type) then that object is
    returned. Otherwise an exception is thrown. */
MarketObjectSP MarketData::GetData(const string&         name,
                                   const CClassConstSP&  type) const{
    return GetObjectData(name, type, 0);
}

MarketObjectSP MarketData::GetObjectData(const string&         name,
                                         const CClassConstSP&  type,
                                         const IModel*         model) const{
    MarketObjectConstSP marketObj;
    try {
        // problems with gcc on solaris - put this statement in a separate
        // try catch block
         marketObj = getData(name, type, model);
    } catch (exception&){
        throw;
    }
    try{
        MarketObjectSP clonedObj = MarketObjectSP(marketObj.clone());
        return clonedObj;
    } catch (exception& e){
        throw ModelException(e, "MarketData::GetObjectData", 
                             "Failed to copy object with "
                             "name "+name+" and type "+
                             marketObj->getClass()->getName());
    }
}

/** Returns the [concrete] type of the market data requested - fails
    if the data does not exist */
CClassConstSP MarketData::getDataType(const string&         name,
                                      const CClassConstSP&  type) const{
    return (getData(name, type)->getClass());
}

/** Returns true if an object of the given name and type exists */
bool MarketData::hasData(const string&         name,
                         const CClassConstSP&  type) const{
    static const string routine("MarketData::hasData");

    bool hasObject = false;

    typedef Imp::CDataMap::const_iterator I;
    // search all entries with the same name to see if we can get a
    // unique match
    MarketObjectConstSP marketObj;
    pair<I,I> limits = my->DataCache.equal_range(name);
    for (I iter = limits.first; !(iter == limits.second) && hasObject == false;
         ++iter){
        if (type->isInstance(iter->second.get())){
            // found a candidate
            hasObject = true;
        }
    }
    return hasObject;
}


/** Asks the supplied model to return a [deep] copy of the
    requested MarketObject object from this cache. The getMarket method
    is then invoked on the returned MarketObject */
MarketObjectSP MarketData::GetData(const IModel*         model,
                                   const string&         name,
                                   const CClassConstSP&  type) const{
    /* ask the model to get the data of the right precise type (eg
       what type of vol exactly) */
    try{
        MarketObjectSP marketObj;
        try{
            // problems with gcc on solaris - put this statement in a separate
            // try catch block
            marketObj = model->GetMarket(this, name, type);
        } catch (exception&){
            throw;
        }
        if (marketObj.get()){ /* occasionally certain market data is not needed
                                 or a hard coded set is going to be used eg
                                 choice of q-smile for cds quanto adjustment.
                                 It seems cleaner to allow the model to return
                                 null in this case */
            // then get the market object to fill in the bits it needs
            model->getComponentMarketData(this, marketObj);
        }
        return marketObj;
    } catch (exception& e){
        throw ModelException(e, "MarketData::GetData", "Failed to get "
                             "market data for model "
                             +model->getClass()->getName());
    }
}

/** Return the reference date (aka value date, today) associated
    with this data */
DateTime MarketData::GetReferenceDate() const{
    return my->RefDate;
}

/** Populates the supplied dateTime with the value date from the
    market data cache. Fails if no date is in the cache and the
    supplied date is empty or if the supplied date is not empty
    and does not match the one in the cache */
void MarketData::GetReferenceDate(DateTime& valueDate) const{
    static const string routine("MarketData::GetReferenceDate");
    if (my->RefDate.empty()){
        if (valueDate.empty()){
            throw ModelException(routine, "No reference date in cache!");
        }

        // set market value date to the first non-empty one
        my->RefDate = valueDate; 
    } else {
        if (valueDate.empty()){
            valueDate = my->RefDate; //return ref date
        } else if (!valueDate.equals(my->RefDate)){
            throw ModelException(routine, "Supplied value date "+
                                 valueDate.toString()+" does not match"
                                 " reference date "+my->RefDate.toString()+
                                 " in cache");
        }
    }
}            

IObject* MarketData::defaultMarketData(){
    return new MarketData();
}


/** Returns the name identifying a correlation given the names of the two
    objects whose correlation is required */
string MarketData::getCorrelationName(const string& id1, 
                                      const string& id2) const{
    Imp::StringPair stringPair(id1, id2);
    string correlationName;
    Imp::SymStringPairHash::const_iterator iter = 
        my->correlationsHash.find(stringPair);
    if (iter == my->correlationsHash.end()){
        throw ModelException("MarketData::getCorrelationName", 
                             "No correlation available with respect to "+
                             id1+" and "+id2);
    }
    return iter->second;
}

/** Returns the name identifying a correlation given the names of the two
    objects whose correlation is required */
bool MarketData::hasCorrelationData(const string& id1, 
                                    const string& id2) const{
    Imp::StringPair stringPair(id1, id2);
    Imp::SymStringPairHash::const_iterator iter = 
        my->correlationsHash.find(stringPair);
    return (iter != my->correlationsHash.end() );
}

/** Sets the name identifying a correlation given the names of the two
    objects identifying the correlation */
void MarketData::setCorrelationName(const string& id1,
                                    const string& id2,
                                    const string& correlationName){
    Imp::StringPair stringPair(id1, id2);
    my->correlationsHash[stringPair] = correlationName;
}

/** Returns the name identifying a correlationTerm given the names of the two
    objects whose correlationTerm is required */
string MarketData::getCorrelationTermName(  const string& id1, 
                                            const string& id2) const{
    Imp::StringPair stringPair(id1, id2);
    string correlationTermName;
    Imp::SymStringPairHash::const_iterator iter = 
        my->correlationTermHash.find(stringPair);
    if (iter == my->correlationTermHash.end()){
        throw ModelException("MarketData::getCorrelationTermName", 
                             "No correlationTerm available with respect to "+
                             id1+" and "+id2);
    }
    return iter->second;
}

/** Returns the name identifying a correlationTerm given the names of the two
    objects whose correlationTerm is required */
bool MarketData::hasCorrelationTermData(const string& id1, 
                                        const string& id2) const{
    Imp::StringPair stringPair(id1, id2);
    Imp::SymStringPairHash::const_iterator iter = 
        my->correlationTermHash.find(stringPair);
    return (iter != my->correlationTermHash.end() );
}

/** Sets the name identifying a correlationTerm given the names of the two
    objects identifying the correlationTerm */
void MarketData::setCorrelationTermName(const string& id1,
                                        const string& id2,
                                        const string& correlationTermName){
    Imp::StringPair stringPair(id1, id2);
    my->correlationTermHash[stringPair] = correlationTermName;
}

/** Returns the name identifying a corrSwapBasisAdj given the names of the two
	regions for which the basis is required */
string MarketData::getCorrSwapSamplingAdjName(const string& id1, 
											  const string& id2) const{
	Imp::StringPair stringPair(id1, id2);	
	Imp::SymStringPairHash::const_iterator iter = 
		my->corrSwapSamplingAdjHash.find(stringPair);
	if (iter == my->corrSwapSamplingAdjHash.end()){
		throw ModelException("MarketData::getCorrSwapSamplingAdjName", 
			"No corrSwapSamplingAdj available with respect to "+
			id1+" and "+id2);
	}
	return iter->second;
}

bool MarketData::hasCorrSwapSamplingAdjData(const string& id1,
											const string& id2) const {
	Imp::StringPair stringPair(id1, id2);
	Imp::SymStringPairHash::const_iterator iter = 
		my->corrSwapSamplingAdjHash.find(stringPair);
	return (iter != my->corrSwapSamplingAdjHash.end() );
												
}

/** Sets the name identifying a corrSwapSamplingAdj given the names of the two
regions for which the sampling adj is required  */
void MarketData::setCorrSwapSamplingAdjName(const string& id1,
											const string& id2,
											const string& corrSwapBasisAdjName) {
	Imp::StringPair stringPair(id1, id2);
	my->corrSwapSamplingAdjHash[stringPair] = corrSwapBasisAdjName;

}

/** Returns the name identifying an fx asset given the names of the two
    yield curves for the two currencies making up the fx rate */
string MarketData::getFXName(const string& riskCcy, 
                             const string& baseCcy) const{
    Imp::StringPair stringPair(riskCcy, baseCcy);
    string fxName;
    Imp::ASymStringPairHash::const_iterator iter = my->fxHash.find(stringPair);
    if (iter == my->fxHash.end()){
        throw ModelException("MarketData::getFXName", 
                             "No fx asset available with respect to risk ccy "+
                             riskCcy+" and base ccy "+baseCcy);
    }
    return iter->second;
}

/** Record what the iso code is for a specific yield curve */
void MarketData::setYieldCurveISOCode(const string& ycName,
                                      const string& ycISOCode){
    my->isoCodesByYC[ycName] = ycISOCode;
}

/** Get what the iso code is for a specific yield curve. This is useful
    if you want to choose what type of data to get inside a yield curve
    dependent upon the iso code */
const string& MarketData::getYieldCurveISOCode(const string& ycName) const{
    Imp::StringTable::const_iterator iter = my->isoCodesByYC.find(ycName);
    if (iter == my->isoCodesByYC.end()){
        throw ModelException("MarketData::getYieldCurveISOCode", "No record "
                             "of yield curve "+ycName+" in cache");
    }
    return iter->second;
}

/** Returns the name identifying a correlation given the names of the two
    objects whose correlation is required */
bool MarketData::hasFXData(const string& riskCcy, 
                           const string& baseCcy) const{
    Imp::StringPair stringPair(riskCcy, baseCcy);
    Imp::ASymStringPairHash::const_iterator iter = my->fxHash.find(stringPair);
    return ( iter != my->fxHash.end() );
}

/** Sets the name identifying a correlation given the names of the two
    objects identifying the correlation */
void MarketData::setFXName(const string& riskCcy,
                           const string& baseCcy,
                           const string& fxName){
    Imp::StringPair stringPair(riskCcy, baseCcy);
    my->fxHash[stringPair] = fxName;
}

/** Returns all ObservableHistorys in the cache for a given
    market observable name and type*/
MarketObjectArraySP MarketData::getAllObservableHistories(const IModel* model,
                                             const string& obsName,
                                             const CClassConstSP&  type) const{
    static const string routine("MarketData::getAllObservableHistories");
    try{
        MarketObjectArraySP marketObjects(new MarketObjectArray(0));
        typedef Imp::ObsHistMap::const_iterator I;
        // search all entries with the same name to see if we can get a match
        pair<I,I> limits = my->obsHistory.equal_range(obsName);
        for (I iter = limits.first; !(iter == limits.second); ++iter){
            // now need to retrieve the observable history
            if (type->getName() == iter->second->histType->getName()) {
                MarketObjectSP obj = GetData(iter->second->cacheName, type);
                obj->getMarket(model, this);
                marketObjects->push_back(obj);
            }
        }
        return marketObjects;
    } catch (exception& e){
        throw ModelException(e, routine, "Failed to retrieve histories for "
                             + obsName);
    }
}

/** Checks whether we already have an ObservableHistory for the given
    Observable name, type and source - note this also throws an Exception
    if we try to add a history of a different type - we assume we can't
    end uf with two histories for different sources for the same observable
    but with the same name in the cache */
bool MarketData::hasObservableHistoryData(const string& obsName,
                                          const string& source,
                                          const CClassConstSP&  type) const{
    static const string routine("MarketData::hasObservableHistoryData");
    try{
        typedef Imp::ObsHistMap::const_iterator I;
        // search all entries with the same name to see if we can get a match
        pair<I,I> limits = my->obsHistory.equal_range(obsName);
        for (I iter = limits.first; !(iter == limits.second); ++iter){
            // first check we haven't already got a different type
            if (!(type->getName() == iter->second->histType->getName())) {
                throw ModelException(routine, "Cannot have two observable "
                            "histories for the same observable ("+obsName+
                            ") but with different types "+type->getName()+
                            " and "+iter->second->histType->getName());
            }
            if (source == iter->second->source) {
                return true;
            }
        }
        return false;
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}

/** Returns the name identifying the ObservableHistory in the cache for a given
    Observable name, type and source - we've already checked we have one
    using hasObservableHistoryData so just get it*/
string MarketData::getObservableHistoryName(const string& obsName,
                                            const string& source,
                                            const CClassConstSP&  type) const{
    static const string routine("MarketData::getObservableHistoryName");
    try{
        typedef Imp::ObsHistMap::const_iterator I;
        // search all entries with the same name to see if we can get a match
        pair<I,I> limits = my->obsHistory.equal_range(obsName);
        for (I iter = limits.first; !(iter == limits.second); ++iter){
            if (source == iter->second->source &&
                type->getName() == iter->second->histType->getName()) {
                return iter->second->cacheName;
            }
        }
        throw ModelException(routine, 
                             "Internal error - no observable history available "
                             "with respect to " +obsName+" and "+source);
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}

/** Sets the name identifying a correlation given the Observable name */
void MarketData::setObservableHistoryName(const string& obsName,
                                          const string& source,
                                          const CClassConstSP&  type,
                                          const string& cacheName){
    try {
        // insert new entry
        // note we allow duplicates for different types and sources
        // and this is called from initialise() which has already
        // checked we haven't seen this one before
        Imp::ObsHistoryDetailsSP details(new Imp::ObsHistoryDetails(source, 
                                                                  cacheName,
                                                                  type));
        Imp::ObsHistMap::value_type entry(obsName, details);
        my->obsHistory.insert(entry);
    } catch (exception& e){
        throw ModelException(e, "MarketData::setObservableHistoryName");
    }
}

/** Fetches any market data required by eg. Business252 */
void MarketData::fetchForNonMarketObjects(IObjectSP     obj, 
                                          IModelConstSP model,
                                          const string& className) const
{
    class Action : public ObjectIteration::IAction
    {
    public:
        // called by ObjectIteration
        bool invoke(ObjectIteration::State& state, IObjectConstSP constObj)
        {
            // Need non-const access to the object which contains market data
            IObjectSP obj = state.getObject();

            IGetMarket& objWithMarketData = (*dynamic_cast<IGetMarket*>(obj.get()));

            // Now get the object to retrieve market data from the cache
            objWithMarketData.getMarket(model, marketData);
            return true;
        }

        const IModel* model;
        const MarketData* marketData;
    };

    // Create the action object
    Action action;
    action.model      = model.get();
    action.marketData = this;
    
    // Build an instance of the class that drives the recursion
    // this is very hacky - and is a temporary fix until we refactor how
    // market data is retrieved.
    ObjectIteration iteration(CClass::forName(className));
        
    iteration.recurse(action, obj);
}

// the market iterator
class MarketIterator: public CObject,
                      public virtual IMap::IIterator {
public:
    static CClassConstSP const TYPE;
    // are there any elements left to iterate over
    virtual bool hasMoreElements() const {      
        return (iter != mkt->my->DataCache.end());
    }

    // get the key for the current element
    virtual const string& getKey() const {
        return iter->first;
    }

    // get the current element (returns a reference)
    virtual IObjectSP getElement() const {
        return IObjectSP(iter->second);
    }
    // increment the iterator
    virtual void increment() {
        ++iter;
    }

    MarketIterator(const MarketDataSP& mkt):
        CObject(TYPE), iter(mkt->my->DataCache.begin()), mkt(mkt) {}
    
    static void load(CClassSP& clazz){
        REGISTER(MarketIterator, clazz);
        SUPERCLASS(CObject)
        IMPLEMENTS(IMap::IIterator);
        EMPTY_SHELL_METHOD(defaultMarketIterator);
    }

    static IObject* defaultMarketIterator(){
        return new MarketIterator();
    }

private:
    MarketIterator():CObject(TYPE) {};
    // data
    MarketData::Imp::CDataMap::const_iterator iter; // $unregistered
    MarketDataSP                              mkt; // $unregistered
};

CClassConstSP const MarketIterator::TYPE = CClass::registerClassLoadMethod(
    "MarketIterator", typeid(MarketIterator), MarketIterator::load);

// Builds an iterator
IMap::IIteratorSP MarketData::createIterator() {    
    return IMap::IIteratorSP(new MarketIterator(MarketDataSP::attachToRef(this)));
}

/** Is this object truly a map ie does toObject()/toMap() return this */
bool MarketData::isTrueMap() const {
    return false;
}

/** Add data to the market */
void MarketData::put(const string& key, const IObjectSP& value) {
    static const string method("MarketData::put");
    try {
        // what if they pass back in an array (see get())?
        // cycle through entries and add them in one at a time?
        if (IArray::TYPE->isInstance(value)){
            IArray* array = DYNAMIC_CAST(IArray, value.get());
            for (int i = 0; i < array->getLength(); i++) {
                put(key, array->get(i));
            }
        } else {
            // must support taking in public representations of objects
            IObjectSP privateRep(CObject::convertToPrivateRep(value));
            MarketObject* mo = dynamic_cast<MarketObject*>(privateRep.get());
            if (!mo) {
                throw ModelException(method, "Object of type "+
                                     privateRep->getClass()->getName()+
                                     " is not a market object");
            }
            // as the key is somewhat redundant here, allow user to skip it
            // if they give it, it has to match though
            if (!key.empty() && key != mo->getRootName()) {
                throw ModelException(method, "market object name ("
                                     + mo->getRootName() + 
                                     ") is inconsistent with key (" + key +")");
            }
            
            AddData(MarketObjectSP(mo));
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Get data from the market - returns an array */
IObjectSP MarketData::get(const string& key) const {
    return GetAllDataWithName(key);
}

// addin support
class MarketDataAddin: public CObject{
    static CClassConstSP const TYPE;

    /* two addin parameters - value date and an
       array of objects to add to the MarketData cache */
    DateTime                 valueDate;
    MarketObjectArray        objects;    // the objects themselves

    /** create a market data cache */
    static IObjectSP create(MarketDataAddin *params){
        static const string routine("MarketDataAddin::create");
        // create market cache with supplied reference date
        CMarketDataSP marketCache(new MarketData(params->valueDate));
        // then add supplied objects to it
        MarketObjectArray& objects = params->objects;
        for (int i = 0; i < objects.size(); i++){
            if (!objects[i]){
                //skip
            } else {
                // use a reference or a copy? A reference seems to make most
                // sense
                marketCache->AddData(objects[i]);
            }
        }
        return marketCache;
    }
 
    MarketDataAddin(): CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(MarketDataAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultMarketDataAddin);
        FIELD(valueDate, "Value date");
        FIELD(objects, 
                     "Array of objects to add to the market data cache");
        FIELD_MAKE_OPTIONAL(objects);
        // register addin for adding objects to cache
        Addin::registerClassObjectMethod(
            "MARKET",
            Addin::MARKET,
            "Creates market data cache",
            TYPE,
            true,
            Addin::returnHandle,
            (Addin::ObjMethod*)create);

    }

    static IObject* defaultMarketDataAddin(){
        return new MarketDataAddin();
    }
};

CClassConstSP const MarketDataAddin::TYPE= CClass::registerClassLoadMethod(
    "MarketDataAddin", typeid(MarketDataAddin), MarketDataAddin::load);

class AddMarketAddin: public CObject{
    static CClassConstSP const TYPE;

    /* two addin parameters - MarketData cache & data to add*/
    CMarketDataSP     market;
    MarketObjectArray objects;    // the objects themselves

    /** clone market data cache and add to it */
    static IObjectSP add(AddMarketAddin *params){
        static const string routine("AddMarketAddin::add");
        // CLONE THE INPUT MARKET
        // Why? A number of reasons:
        // (a) that's what the method says that it does
        // (b) it's used as an Excel addin, and cloning and returning a new
        // market handle is the only sensible behaviour otherwise you are
        // changing the original without updating the handle, so you get 
        // horrible calculation order dependent behaviour.
        CMarketDataSP market(copy(params->market.get()));

        MarketObjectArray& objects = params->objects;
        for (int i = 0; i < objects.size(); i++){
            if (!objects[i]){
                //skip
            } else {
                // use a reference or a copy? A reference seems to make most
                // sense
                market->AddData(objects[i]);
            }
        }
        return market; // this is a cloned and extended version of the original
    }
 
    AddMarketAddin(): CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(AddMarketAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultAddMarketAddin);
        FIELD(market, "market data cache");
        FIELD(objects, 
                     "Array of objects to add to the market data cache");
        // register addin for adding objects to cache
        Addin::registerClassObjectMethod(
            "EXTEND_MARKET",
            Addin::MARKET,
            "copy market data cache and add to it",
            TYPE,
            true,
            Addin::returnHandle,
            (Addin::ObjMethod*)add);

    }

    static IObject* defaultAddMarketAddin(){
        return new AddMarketAddin();
    }
};

CClassConstSP const AddMarketAddin::TYPE= CClass::registerClassLoadMethod(
    "AddMarketAddin", typeid(AddMarketAddin), AddMarketAddin::load);

// addin support
class GetMarketDataAddin: public CObject{
    static CClassConstSP const TYPE;

    CMarketDataSP        object;

    /** create a market data cache */
    static IObjectSP getObjects(GetMarketDataAddin* params){
        static const string routine("GetMarketDataAddin::getObjects");

        typedef multimap<const string, MarketObjectSP> CDataMap;

        smartPtr<MarketObjectArray> marketObjArray(new MarketObjectArray(0));
        for (CDataMap::const_iterator myIterator = 
                 params->object->my->DataCache.begin();
             !(myIterator == params->object->my->DataCache.end());
             ++myIterator){
                 marketObjArray->push_back(myIterator->second);
                 // do nothing for the time being
        }
        return marketObjArray;
    }
 
    GetMarketDataAddin(): CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(GetMarketDataAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultGetMarketDataAddin);
        FIELD(object, "Handle to the Market Data Cache");
        // register addin for adding objects to cache
        Addin::registerClassObjectMethod(
            "GET_MARKET",
            Addin::MARKET,
            "Retrieves all objects from the Market Data Cache",
            TYPE,
            true,
            Addin::returnHandle,
            (Addin::ObjMethod*)getObjects);

    }

    static IObject* defaultGetMarketDataAddin(){
        return new GetMarketDataAddin();
    }
};

CClassConstSP const GetMarketDataAddin::TYPE= CClass::registerClassLoadMethod(
    "GetMarketDataAddin", typeid(GetMarketDataAddin), GetMarketDataAddin::load);

// addin support
class GetReferenceDateAddin: public CObject {
    static CClassConstSP const TYPE;

    CMarketDataSP        object;

    /** create a market data cache */
    static IObjectSP getRefDate(GetReferenceDateAddin* params){
        static const string routine("GetReferenceDateAddin::getRefDate");

        return IObjectSP(params->object->GetReferenceDate().clone());
    }
 
    GetReferenceDateAddin(): CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(GetReferenceDateAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultReferenceDateAddin);
        FIELD(object, "Handle to the Market Data Cache");
        // register addin for adding objects to cache
        Addin::registerClassObjectMethod(
            "GET_REFERENCE_DATE",
            Addin::MARKET,
            "returns a handle to the reference date in the Market Data Cache",
            TYPE,
            true,
            Addin::returnHandle,
            (Addin::ObjMethod*)getRefDate);

    }

    static IObject* defaultReferenceDateAddin(){
        return new GetReferenceDateAddin();
    }
};

CClassConstSP const GetReferenceDateAddin::TYPE= 
CClass::registerClassLoadMethod(
    "GetReferenceDateAddin", typeid(GetReferenceDateAddin), 
    GetReferenceDateAddin::load);

// addin support for merge markets together
typedef array<MarketDataSP, MarketData> MarketDataArray;
DEFINE_TEMPLATE_TYPE(MarketDataArray);

typedef smartPtr<MarketDataArray> MarketDataArraySP;
class MergeMarketsAddin: public CObject{
    static CClassConstSP const TYPE;

    MarketDataArraySP   markets;

    /** create a market data cache */
    static IObjectSP mergeMarkets(MergeMarketsAddin* params){
        return params->run();
    }
    IObjectSP run(){
        static const string routine("MergeMarketsAddin::run");
        MarketDataSP all;
        if (markets->empty()){
            all = MarketDataSP(new MarketData());
        }
        for (int i = 0; i < markets->size(); i++){
            if (!(*markets)[i]){
                throw ModelException(routine, "Market number "+
                                     Format::toString(i+1)+" is null");
            }
            const DateTime& refDate = (*markets)[i]->GetReferenceDate();
            if (i == 0){
                all = MarketDataSP(new MarketData(refDate));
            } else if (!refDate.empty()){
                const DateTime& curRefDate = all->GetReferenceDate();
                if (!refDate.equals(curRefDate)){
                    throw ModelException(routine, 
                                         "All markets must have the same"
                                         " reference date");
                }
            }
            // then merge current market into 'all'
            all->addData(*(*markets)[i]);
        }
        return all;
    }
 
    MergeMarketsAddin(): CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(MergeMarketsAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultMergeMarketsAddin);
        FIELD(markets, "Array of markets to merge");
        // register addin for adding objects to cache
        Addin::registerClassObjectMethod(
            "MERGE_MARKETS",
            Addin::MARKET,
            "Merges array of markets into one big market",
            TYPE,
            true,
            Addin::returnHandle,
            (Addin::ObjMethod*)mergeMarkets);
    }

    static IObject* defaultMergeMarketsAddin(){
        return new MergeMarketsAddin();
    }
};

CClassConstSP const MergeMarketsAddin::TYPE= CClass::registerClassLoadMethod(
    "MergeMarketsAddin", typeid(MergeMarketsAddin), MergeMarketsAddin::load);


// a proxy interface for MarketData so fits in with Global DR Interface
// for serialisation etc.

// this represents the cache portion of the proxy interface
class MarketCache: public CObject,
                   virtual public IReadableMap,
                   virtual public IWriteableMap {
public:
    static CClassConstSP const TYPE;

    // Map interfaces
    // Builds an iterator
    virtual IMap::IIteratorSP createIterator() {
        return mkt->createIterator();
    }

    /** Is this object truly a map ie does toObject()/toMap() return this */
    virtual bool isTrueMap() const {
        return true;
    }

    /** Add data to the market */
    virtual void put(const string& key, const IObjectSP& value) {
        // what if they pass back in an array (see get())?
        // cycle through entries and add them in one at a time?
        IArray* array = dynamic_cast<IArray*>(value.get());
        if (!array) {
            mkt->put(key, value);
        }
        else {
            for (int i = 0; i < array->getLength(); i++) {
                mkt->put(key, array->get(i));
            }
        }
    }

    /** Get data from the market - returns an array */
    virtual IObjectSP get(const string& key) const {
        return mkt->get(key);
    }

    // over-ride clone - default would lose the market as not registered
    IObject* clone() const {
        MarketCache* clone = new MarketCache(mkt.get());
        return clone;
    }

    MarketCache(MarketData* market): CObject(TYPE), mkt(market) {}

private:
    MarketDataSP mkt;  // I'm really a market in disguise $unregistered

    friend class MarketDataProxy;  // oh I'm so lazy
    /** for reflection */
    MarketCache(): CObject(TYPE), mkt(new MarketData) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); 
        REGISTER(MarketCache, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IReadableMap);
        IMPLEMENTS(IWriteableMap);
        EMPTY_SHELL_METHOD(defaultMarketCache);
    }

    static IObject* defaultMarketCache(){
        return new MarketCache();
    }
};

CClassConstSP const MarketCache::TYPE = CClass::registerClassLoadMethod(
    "MarketCache", typeid(MarketCache), MarketCache::load);

typedef smartPtr<MarketCache> MarketCacheSP;

class MarketDataProxy: public CObject {
public:
    static CClassConstSP const TYPE;

private:
    DateTime      today;
    MarketCacheSP cache;

    typedef smartPtr<MarketDataProxy> MarketDataProxySP;

    static IObjectSP toMarketData(const IObjectConstSP& obj) {
        const MarketDataProxy* proxy = dynamic_cast<const MarketDataProxy*>(obj.get());
        if (!proxy) {
            throw ModelException("MarketDataProxy::toMarketData",
                                 "object is not a MarketDataProxy");
        }
    
        // if there's something in the cache, then that's the market
        MarketDataSP market;

        if (!proxy->cache.get()) {
            market = MarketDataSP(new MarketData(proxy->today));
        }
        else {
            // excuse the cast, but it really is a MarketData
            MarketData* copy = dynamic_cast<MarketData*>(proxy->cache->mkt->clone());
            market = MarketDataSP(copy);
            if (market->GetReferenceDate().empty()) {
                // market was set up without a reference date
                market->my->RefDate = proxy->today;
            }
        }
        return market;
    }

    static IObjectSP fromMarketData(const IObjectConstSP& obj) {
        const MarketData* market = dynamic_cast<const MarketData*>(obj.get());
        if (!market) {
            throw ModelException("MarketDataProxy::fromMarketData",
                                 "object is not a MarketData");
        }
    
        MarketData* mkt = const_cast<MarketData*>(market);

        MarketDataProxySP proxy(new MarketDataProxy(mkt->GetReferenceDate()));
        // stuff ourself in as the proxy's cache
        proxy->cache = MarketCacheSP(new MarketCache(mkt));
        return proxy;
    }

    MarketDataProxy(const DateTime& today) : CObject(TYPE),today(today) {}

    /** for reflection */
    MarketDataProxy():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); 
        REGISTER(MarketDataProxy, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultMarketDataProxy);
        registerObjectProxy(MarketData::TYPE,
                            MarketDataProxy::TYPE,
                            fromMarketData,
                            toMarketData);
        FIELD(today, "today");
        FIELD(cache, "cache");
        FIELD_MAKE_OPTIONAL(cache);
    }

    static IObject* defaultMarketDataProxy(){
        return new MarketDataProxy();
    }  
};

CClassConstSP const MarketDataProxy::TYPE = CClass::registerClassLoadMethod(
    "MarketDataProxy", typeid(MarketDataProxy), load);

void MarketData::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    clazz->setDRIProxyType(MarketDataProxy::TYPE); // use proxy for dri
    REGISTER(MarketData, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IReadableMap);
    IMPLEMENTS(IWriteableMap);
    EMPTY_SHELL_METHOD(defaultMarketData);
}

CClassConstSP const MarketData::TYPE = CClass::registerClassLoadMethod(
    "MarketData", typeid(MarketData), load);

/**
 * Test for MarketData::plugin() mechanism
 */

struct MarketDataPluginCacheTest: CObject, virtual ClientRunnable {

    static const CClassConstSP TYPE;
    static IObject* emptyShell();
    static void load(CClassSP&);

    struct Thing: MarketObject {

        static const CClassConstSP TYPE;
        static IObject* emptyShell();
        static void load(CClassSP&);

        string name; // $required

        Thing(const string& name = ""):
            MarketObject(TYPE),
            name(name)
        {}

        static MarketObjectSP SP(const string& name) {
            return MarketObjectSP(new Thing(name));
        }

        string getName() const {
            return name;
        }
    };

    struct ThingNames: CObject {

        static const CClassConstSP TYPE;
        static IObject* emptyShell();
        static void load(CClassSP&);

        string str; // $required

        static int count;

        ThingNames(): CObject(TYPE) {}

        ThingNames(const MarketData& market):
            CObject(TYPE)
        {
            MarketObjectArraySP them = market.GetAllDataWithType(Thing::TYPE);
            for (int t = 0; t < them->size(); ++t) {
                str += (*them)[t]->getName() + " ";
            }

            str += "#" + Format::toString(count);
            ++count;
        }

        static IObjectConstSP SP(const MarketData& market) {
            return IObjectConstSP(new ThingNames(market));
        }

        string toString() const {
            return str;
        }
    };

    MarketDataPluginCacheTest(): CObject(TYPE) {}

    IObjectSP run() {
        try {
            ThingNames::count = 0;

            ostringstream o;

#           define DO(X) o << "  " << #X << ";\n"; X;
#           define CHECK(X, V)                          \
              do {                                      \
                string v((X).toString());               \
                o << "  " << #X << " = " << v;           \
                if (v != (V)) o << " (THAT'S A BUG)";   \
                o << "\n";                              \
              }                                         \
              while (0)
#           define COMMENT(C) o << "\n" << (C) << "\n\n";

            o << "----\n\n";

            {
                MarketDataSP market(new MarketData());
                DO( market->AddData(Thing::SP("Alice")) );
                DO( market->AddData(Thing::SP("Bob")) );

                COMMENT( "Trigger creation of a ThingNames cache:" );
                CHECK( *market->plugin<ThingNames>(), "Alice Bob #0" );

                COMMENT( "A cloned MarketData shares the same cache (#0):" );
                DO( MarketDataSP market2(dynamic_cast<MarketData*>(market->clone())) );
                CHECK( *market2->plugin<ThingNames>(), "Alice Bob #0" );

                COMMENT( "Altering the original causes it make a new cache:" );
                DO( market->AddData(Thing::SP("Chloe")) );
                CHECK( *market->plugin<ThingNames>(), "Alice Bob Chloe #1" );
                CHECK( *market2->plugin<ThingNames>(), "Alice Bob #0" );
            }
            
            o << "\n----\n\n";

            {
                MarketDataSP market(new MarketData());
                DO( market->AddData(Thing::SP("Alice")) );
                DO( market->AddData(Thing::SP("Bob")) );

                COMMENT( "The cache is shared even if we clone before asking for it:" );
                DO( MarketDataSP market2(dynamic_cast<MarketData*>(market->clone())) );
                CHECK( *market->plugin<ThingNames>(), "Alice Bob #2" );
                CHECK( *market2->plugin<ThingNames>(), "Alice Bob #2" );
            }

            return IObjectSP(CString::create(o.str()));
        }
        catch (exception& e) {
            throw ModelException(e, __FUNCTION__);
        }
    }
};

int MarketDataPluginCacheTest::ThingNames::count = 0;

IObject* MarketDataPluginCacheTest::ThingNames::emptyShell() {
  return new MarketDataPluginCacheTest::ThingNames();
}

void MarketDataPluginCacheTest::ThingNames::load(CClassSP& clazz) {
  REGISTER(MarketDataPluginCacheTest::ThingNames, clazz);
  SUPERCLASS(CObject);
  FIELD(str, "str");
  EMPTY_SHELL_METHOD(emptyShell);
}

typedef MarketDataPluginCacheTest::ThingNames MarketDataPluginCacheTest_ThingNames;
CClassConstSP const MarketDataPluginCacheTest_ThingNames::TYPE = CClass::registerClassLoadMethod(
  "MarketDataPluginCacheTest::ThingNames", typeid(MarketDataPluginCacheTest_ThingNames), MarketDataPluginCacheTest_ThingNames::load);

IObject* MarketDataPluginCacheTest::Thing::emptyShell() {
  return new MarketDataPluginCacheTest::Thing();
}

void MarketDataPluginCacheTest::Thing::load(CClassSP& clazz) {
  REGISTER(MarketDataPluginCacheTest::Thing, clazz);
  SUPERCLASS(MarketObject);
  FIELD(name, "name");
  EMPTY_SHELL_METHOD(emptyShell);
}

typedef MarketDataPluginCacheTest::Thing MarketDataPluginCacheTest_Thing;
CClassConstSP const MarketDataPluginCacheTest_Thing::TYPE = CClass::registerClassLoadMethod(
  "MarketDataPluginCacheTest::Thing", typeid(MarketDataPluginCacheTest_Thing), MarketDataPluginCacheTest_Thing::load);

CClassConstSP const MarketDataPluginCacheTest::TYPE = CClass::registerClassLoadMethod(
  "MarketDataPluginCacheTest", typeid(MarketDataPluginCacheTest), MarketDataPluginCacheTest::load);

IObject* MarketDataPluginCacheTest::emptyShell() {
  return new MarketDataPluginCacheTest();
}

void MarketDataPluginCacheTest::load(CClassSP& clazz) {
  REGISTER(MarketDataPluginCacheTest, clazz);
  SUPERCLASS(CObject);
  EMPTY_SHELL_METHOD(emptyShell);
}

DRLIB_END_NAMESPACE
