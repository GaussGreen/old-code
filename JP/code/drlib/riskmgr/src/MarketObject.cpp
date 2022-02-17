//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MarketObject.cpp
//
//   Description : Base class for objects stored in market data cache
//
//   Author      : Mark A Robson
//
//   Date        : 19 Mar 2001
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#define QLIB_MARKETOBJECT_CPP
#include "edginc/MarketObject.hpp"
#include "edginc/IModel.hpp"
#include "edginc/MarketData.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/Writer.hpp"
#include "edginc/IMarketObjectQualifier.hpp"

DRLIB_BEGIN_NAMESPACE

INamedObject::INamedObject() {}
INamedObject::~INamedObject() {}

void INamedObject::load(CClassSP& clazz){
    REGISTER_INTERFACE(INamedObject, clazz);
    EXTENDS(IObject);
}

CClassConstSP const INamedObject::TYPE = CClass::registerInterfaceLoadMethod(
    "INamedObject", typeid(INamedObject), load);

MarketObject::~MarketObject(){}

/** Populates the object with the market data that this object
    needs.  This method is invoked as the getMarket chains down
    from the instrument to the specific instance of market
    data. The default implementation provided by MarketObject is
    to do nothing. Market data objects that require other pieces
    of market data (eg an XCB requires assets) need to override
    this method */
void MarketObject::getMarket(const IModel* model, const MarketData* market){
}

/** Returns the name by which this piece of market data is referenced by
    other market data or by instruments/models etc. For example a
    collection of different type of parameterized equity volatilities for
    the FTSE would all have the same name eg "FTSE". Similarly, base and
    swaption interest rate volatilities for the same yield curve would
    have the same root name. The default implementation is to return getName() */
string MarketObject::getRootName() const
{
    // Default implementation for backwards compatibility
    return getName();
}

/** When more than one market object of appropriate type has the same root
    name an additional qualifier is needed in order to determine which
    instance should be used. It is proposed that an object of type
    IMarketObjectQualifier is used to qualify MarketObjects which have the
    same "root name". This type would contain only one method, namely a
    pure virtual equals(IMarketObjectQualifier*) method. One obvious
    implementation of such an interface is a class that just wraps a single
    string. Following the examples in getRootName() this could be, for
    example, "FLOW" or "EXOTIC" for equity volatilies or "BASE" or
    "SWAPTION" for interest rate volatilities. (These strings are only
    examples - there are no predefined values.). The default implementation
    returns a NULL smart pointer. */
IMarketObjectQualifierSP MarketObject::getMarketObjectQualifier() const
{
    // Default MarketObject has no IMarketObjectQualifier
    return IMarketObjectQualifierSP();
}

/** Initialises this piece of market data. It is invoked ONCE only
    - immediately after this object is placed in the cache. The
    default implementation provided by MarketObject is to do
    nothing. */
void MarketObject::initialise(MarketData* market){
    // do nothing
}

/** Converts this object to an instance of the requiredType. Throws an
    exception if a conversion to the required Type is not
    supported. This method supports building the wrapper equivalent of
    a given market object */
void MarketObject::convert(IObjectSP&    object,
                           CClassConstSP requiredType) const{
    if (!MarketObjectWrapper::TYPE->isAssignableFrom(requiredType)){
        throw ModelException("MarketObject::convert", "Can't convert a "
                             "MarketObject into a "+requiredType->getName());
    }
    // create new instance
    IObjectSP newObject(requiredType->newInstance());
    MarketObjectWrapper& wrapper =
        dynamic_cast<MarketObjectWrapper&>(*newObject);
    /* as we have the original smart pointer we can avoid copying the
       market object */
    wrapper.setObject(MarketObjectSP::dynamicCast(object));
    wrapper.useCache = false;
    object = newObject; // return the wrapper as our converted object
}

/** If another object acts as a proxy for this object then return the
    type of the proxy. Eg the VolPreferred object returns the type of
    the vol that it wants to choose. This method is used when pulling
    objects out of the market data cache so that the proxy (eg param vol)
    is used directly instead of the original (eg VolPreferred). */
CClassConstSP MarketObject::proxyType(const MarketData* market,
                                      const IModel*     model) const{
    return 0; // no proxy
}

string MarketObject::toString() const {
    return getClass()->getName() + " \"" + getName() + "\"";
}

MarketObject::MarketObject(const CClassConstSP& clazz): CObject(clazz){
    // empty
}

/** Invoked when Class is 'loaded' */
void MarketObject::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(MarketObject, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(ITypeConvert);
    IMPLEMENTS(IGetMarket);
}

CClassConstSP const MarketObject::TYPE = CClass::registerClassLoadMethod(
    "MarketObject", typeid(MarketObject), MarketObject::load);
DEFINE_TEMPLATE_TYPE(MarketObjectArray);

/// MarketObjectWrapper class
MarketObjectWrapper::~MarketObjectWrapper(){}

MarketObjectWrapper::MarketObjectWrapper(const CClassConstSP& clazz):
    CObject(clazz), useCache(true) {}

MarketObjectWrapper::MarketObjectWrapper(const CClassConstSP& clazz,
                                         const string&        name):
    CObject(clazz), name(name), useCache(true) {}
    
MarketObjectWrapper::MarketObjectWrapper(const CClassConstSP&  clazz, 
                                         bool                  useCache):
    CObject(clazz), useCache(useCache) {}

MarketObjectWrapper::MarketObjectWrapper(const MarketObjectWrapper& rhs):
    CObject(rhs.getClass()), name(rhs.name), useCache(rhs.useCache){}

MarketObjectWrapper& MarketObjectWrapper::operator=(
    const MarketObjectWrapper& rhs){
    name = rhs.name;
    useCache = rhs.useCache;
    return *this;
}

void MarketObjectWrapper::validatePop2Object(){
    /* calculate derived value of useCache (could consider ditching
       useCache field and driving off name.empty() directly) */
    if (name.empty() || getMO().get()){
        useCache = false;
    } else {
        useCache = true;
    }
}

/** Calls CObject::clone and then sets useCache flag */
IObject* MarketObjectWrapper::clone() const{
    MarketObjectWrapper& copy =
        dynamic_cast<MarketObjectWrapper&>(*CObject::clone());
    copy.useCache = useCache;
    return &copy;
}

/** sets the name field within the wrapper. Only valid for uninitialised
    MarketObjectWrappers */
void MarketObjectWrapper::setName(const string& name){
    if (!this->name.empty() || !useCache){
        throw ModelException("MarketObjectWrapper::setName", "Wrapper has "
                             "already been initialised");
    }
    this->name = name;
}

/** Substitutes a name for another, and checks that the object is null. Note 
 * that market data needs to be fetched *after* this substitution for it to have any effect.
 * The return value indicates whether the substitution was successful.
 * An exception is thrown if an object already exists in the wrapper, since
 * we don't know what the correct action is in that case.
 */
bool MarketObjectWrapper::substituteName(const string& oldName, const string& newName){
    if (!!getMO()) {
        throw ModelException("MarketObjectWrapper::substituteName",
                             "Attempt to change the name of a market object which is not null. This is not supported");
    }
    if (this->name == oldName) {
        this->name = newName;
        return true;
    }
    return false;
}

/** Retrieves market data from cache if needed. Asks for data of type of
    supplied clazz to select market data */
void MarketObjectWrapper::getData(const IModel*           model, 
                                  const CMarketDataSP&    market,
                                  const CClassConstSP&    clazz,
                                  const string&           domesticYCName){
    getData(model, market.get(), clazz, domesticYCName);
}

/** Same as above, but no change in domesticYCName */
void MarketObjectWrapper::getData(const IModel*                  model, 
                                  const smartPtr<MarketData>&    market,
                                  const CClassConstSP&           clazz){
    getData(model, market.get(), clazz, "");
}
    

/** Same as above, but no change in domesticYCName */
void MarketObjectWrapper::getData(const IModel*                  model, 
                                  const MarketData*              market,
                                  const CClassConstSP&           clazz){
    getData(model, market, clazz, "");
}

/** Retrieves market data from cache if needed. Asks for data of type of
    supplied clazz to select market data */
void MarketObjectWrapper::getData(const IModel*           model, 
                                  const MarketData*       market,
                                  const CClassConstSP&    clazz,
                                  const string&           domesticYCName){
    MarketObjectSP mo;
    if (useCache){
        if (name.empty()){
            // this might be due to some market data that is optional but
            // is being defaulted
            mo = getMO();
            if (!mo){
                // no name and not defaulted, so fail
                throw ModelException("MarketObjectWrapper::getData", "Name for"
                                     " object is empty when looking for object"
                                     " of type "+clazz->getName());
            }
            // otherwise assume that the data has been defaulted and all is ok
        } else {
            if (domesticYCName.empty()){
                mo = market->GetData(model, name, clazz);
            } else {
                mo = model->getMarketObject(market, name, 
                                            clazz, domesticYCName);
            }
        }
    } else {
        // ensure market object has all its market data
        mo = getMO();
        if (!mo.get())
        {
            throw ModelException("MarketDataFetcher::getComponentMarketData",
                                "Could not get data for MarketObject " + getName() + "!");
        }
        if (domesticYCName.empty()){
            model->getComponentMarketData(market, mo);
        } else {
            model->getComponentMarketData(market, mo, domesticYCName);
        }
    }
    mo = model->modifyMarketData(market, clazz, mo);
    setObject(mo);
}

/** Used for building instances of market object wrappers from
    strings. This function is only called by CString when requiredType
    is a MarketObject (or derived from a MarketObject) */
IObject* MarketObjectWrapper::createFromString(CClassConstSP requiredType,
                                               const string& data)
{
    if (!MarketObjectWrapper::TYPE->isAssignableFrom(requiredType)){
        throw ModelException("MarketObject::convert", requiredType->getName()+
                             " not derived from MarketObject");
    }
    MarketObjectWrapper& wrapper = 
        dynamic_cast<MarketObjectWrapper&>(*requiredType->newInstance());
    wrapper.name = data;
    return &wrapper;
}

/** Returns true if the wrapper is getting its market data from the
    cache */
bool MarketObjectWrapper::usingCache() const{
    return useCache;
}

/** sets whether the wrapper should use the cache to retrieve the corresponding
    market data or not */
void MarketObjectWrapper::setCacheUse(bool useCache){
    this->useCache = useCache;
}

/** Returns the name of the market object that is contained in this
    wrapper. More precisely, if the wrapper is using the cache to get
    the data the original string supplied when the wrapper was created
    is returned otherwise the getName method is called on the object */
string MarketObjectWrapper::getName() const{
    if (useCache){
        return name;
    }
    MarketObjectSP mo(getMO());
    if (!mo){ // avoid crash/going in circles if object is null
        return NO_NAME;
    }
    return mo->getName();
}
  
// is this an empty wrapper i.e. no name and no object
bool MarketObjectWrapper::isEmpty() const {
    return (name.empty() && !getMO().get());
}

//// Routes through equalTo. Method added to support
//// instantiating array template
bool MarketObjectWrapper::operator==(const MarketObjectWrapper& rhs) const{
    return equalTo(&rhs);
}

///////////////////////////////////////////////////////////////////////////////
/** We need it to sort the array of MarketObjects, e.g. MarketWrapper<X>s *////
///////////////////////////////////////////////////////////////////////////////
bool MarketObjectWrapper::operator<(const MarketObjectWrapper & g2Compare)const
{
    const static string routine = "MarketObjectWrapper::operator<";
    try
    {   
        return (this->getName() < g2Compare.getName());
    }
    catch (exception& e)
    {
        throw ModelException(e, routine);
    }
}

/** Either just writes out name (if using cache) or object */
void  MarketObjectWrapper::write(const string& tag, Writer* writer) const{
    if (useCache){
        CStringSP nameObj(CString::create(name));
        nameObj->write(tag, writer);
    } else {
        MarketObjectSP mo(getMO());
        if (!mo){
            writer->writeNull(tag);
        } else {
            mo->write(tag, writer);
        }
    }
}

/** Converts this object to an instance of the
    requiredType. Throws an exception if a conversion to the
    required Type is not supported. This method supports building
    the MarketWrapper<X> from MarketWrapper<Y> when Y is derived from X */
void MarketObjectWrapper::convert(IObjectSP&    object,
                                  CClassConstSP requiredType) const{
    static const string method("MarketObjectWrapper::convert");
    // idiot check
    if (!MarketObjectWrapper::TYPE->isAssignableFrom(requiredType)){
        throw ModelException(method, "Can't convert a "+
                             getClass()->getName()+" into a "+
                             requiredType->getName());
    }
    MarketObjectWrapper& wrapper = 
        dynamic_cast<MarketObjectWrapper&>(*requiredType->newInstance());
    object = IObjectSP(&wrapper); // return the wrapper as our converted object
    // check the wrapped MarketObjects are correct
    if (!wrapper.getMOType()->isAssignableFrom(getMOType())){
        throw ModelException(method, "Can't convert a MarketWrapper for "+
                             getMOType()->getName()+" into a MarketWrapper for "
                             + wrapper.getMOType()->getName());
    }
    // just need to copy over the fields
    wrapper.name = name;
    wrapper.useCache = useCache;
    MarketObjectSP mo(getMO());
    if (mo.get()){
        wrapper.setObject(mo);
    }
}

/** Throws an exception with message getClass()->getName()+ " with name "+
    getName()+" is null")*/
void MarketObjectWrapper::throwExceptionForNullObject() const{
    throw ModelException("MarketObjectWrapper::operator->", 
                         getClass()->getName()+ " with name "+
                         getName()+" is null");
}

/** Invoked when this class is 'loaded' */
void MarketObjectWrapper::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(MarketObjectWrapper, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(ITypeConvert);
    FIELD(name, "Name for piece of market data");
    FIELD_MAKE_OPTIONAL(name);
    CString::registerObjectFromStringMethod(TYPE, createFromString);
}

CClassConstSP const MarketObjectWrapper::TYPE = 
CClass::registerClassLoadMethod("MarketObjectWrapper",
                                typeid(MarketObjectWrapper), load);

// This is the name returned by the getName method if not using the cache
// and there is no embedded object
const string MarketObjectWrapper::NO_NAME = "null";

DRLIB_END_NAMESPACE
