//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MarketObject.hpp
//
//   Description : Base class for objects stored in market data cache
//
//   Author      : Mark A Robson
//
//   Date        : 19 Mar 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_MARKET_OBJECT_HPP
#define EDG_MARKET_OBJECT_HPP

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/Array.hpp"
#include "edginc/TypeConvert.hpp"
#include "edginc/GetMarket.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE
class IModel;
class MarketData;
class MarketObject;
FORWARD_DECLARE(IMarketObjectQualifier)

FORWARD_DECLARE(INamedObject)

/**
 * Interface for objects that have a name
 *
 * This is really for use where you want to use MarketObject::getName() for
 * identification purposes but don't want to assume all the other methods.
 */

class RISKMGR_DLL INamedObject: virtual public IObject {

    static void load(CClassSP& clazz);

public:

    static CClassConstSP const TYPE;

    /** Returns the name of this object. This is the name with which
        it is stored in the market data cache and is the name with
        which results (eg tweaks) should be reported against */
    virtual string getName() const = 0;

    INamedObject();
    ~INamedObject();
};

class RISKMGR_DLL MarketObject: public CObject,
                                virtual public INamedObject,
                                virtual public IGetMarket,
                                virtual public ITypeConvert{
public:
    static CClassConstSP const TYPE;

    virtual ~MarketObject();

    /** Populates the object with the market data that this object
        needs.  This method is invoked as the getMarket chains down
        from the instrument to the specific instance of market
        data. The default implementation provided by MarketObject is
        to do nothing. Market data objects that require other pieces
        of market data (eg an XCB requires assets) need to override
        this method */
    virtual void getMarket(const IModel* model, const MarketData* market);

    /** Returns the name by which this piece of market data is referenced by
        other market data or by instruments/models etc. For example a
        collection of different type of parameterized equity volatilities for
        the FTSE would all have the same name eg "FTSE". Similarly, base and
        swaption interest rate volatilities for the same yield curve would
        have the same root name. The default implementation is to return getName() */
    virtual string getRootName() const;

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
    virtual IMarketObjectQualifierSP getMarketObjectQualifier() const;

    /** Returns the identifier for this object as supplied by the client. It
        is the name typically used to identify sensitivities. Note that it is
        not the name used when retrieving market data (see getRootName()). As
        such it is not a ‘magic’ string, i.e., its value does not have to
        appear anywhere else in the input data. */
    virtual string getName() const = 0;

    /** Initialises this piece of market data. It is invoked ONCE only
        - immediately after this object is placed in the cache. The
        default implementation provided by MarketObject is to do
        nothing. */
    virtual void initialise(MarketData* market);

    /** Converts this object to an instance of the
        requiredType. Throws an exception if a conversion to the
        required Type is not supported. This method supports building
        the wrapper equivalent of a given market object */
    virtual void convert(IObjectSP&    object,
                         CClassConstSP requiredType) const;

    /** If another object acts as a proxy for this object then return
        the type of the proxy (else return null). Eg the VolPreferred
        object returns the type of the vol that it wants to
        choose. This method is used when pulling objects out of the
        market data cache so that the proxy (eg param vol) is used
        directly instead of the original (eg VolPreferred).  The
        default implementation null */
    virtual CClassConstSP proxyType(const MarketData* market,
                                    const IModel*     model) const;

    /** A very brief English name or description of the object,
        for use in constructing e.g. error messages.  This implementation
        returns '<type> "<name>"' */
    virtual string toString() const;

protected:
    MarketObject(const CClassConstSP& clazz);
private:
    static void load(CClassSP& clazz);
    MarketObject(const MarketObject& rhs);
    MarketObject& operator=(const MarketObject& rhs);
};

#ifndef QLIB_MARKETOBJECT_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartConstPtr<MarketObject>);
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartPtr<MarketObject>);
#else
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartConstPtr<MarketObject>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartPtr<MarketObject>);
#endif

// typedef for smart pointers to MarketObject
typedef smartConstPtr<MarketObject> MarketObjectConstSP;
typedef smartPtr<MarketObject> MarketObjectSP;

// support for arrays of MarketObject's
typedef array<MarketObjectSP, MarketObject> MarketObjectArray;
typedef smartPtr<MarketObjectArray> MarketObjectArraySP;
typedef smartConstPtr<MarketObjectArray> MarketObjectArrayConstSP;

class RISKMGR_DLL MarketObjectWrapper: public CObject,
                           public virtual ITypeConvert {
public:
    friend class MarketObject;
    static CClassConstSP const TYPE;
    static const string NO_NAME;

    virtual ~MarketObjectWrapper();

    MarketObjectWrapper(const MarketObjectWrapper& rhs);
    MarketObjectWrapper& operator=(const MarketObjectWrapper& rhs);

    /** Called after an instance of a class (which implements this
        interface) is constructed via 'pop2Object' ie reflection is
        used to fill the internal fields of the object directly.
        Here */
    virtual void validatePop2Object();

    /** Calls CObject::clone and then sets useCache flag */
    virtual IObject* clone() const;

    /** Either just writes out name (if using cache) or object */
    virtual void write(const string& tag, Writer* writer) const;

    /** Converts this object to an instance of the
        requiredType. Throws an exception if a conversion to the
        required Type is not supported. This method supports building
        the MarketWrapper<X> from MarketWrapper<Y> when Y is derived from X */
    virtual void convert(IObjectSP&    object,
                         CClassConstSP requiredType) const;

    /** Returns the Market Object within the wrapper. Does NOT throw
        an exception if null */
    virtual MarketObjectSP getMO() const = 0;

    /** Returns the type of the market object that this is a wrapper for */
    virtual CClassConstSP getMOType() const = 0;

    /** sets the name field within the wrapper. Only valid for uninitialised
        MarketObjectWrappers */
    void setName(const string& name);

    /** Substitutes a name for another, and checks that the object is null. Note 
        that market data needs to be fetched *after* this substitution for it to have any effect.*/
    bool substituteName(const string& oldName, const string& newName);

    /** Stores the supplied object in the wrapper - no copy should be made */
    virtual void setObject(const MarketObjectSP& marketObject) = 0;

    /** Retrieves market data from cache if needed. Asks for data of type of
        supplied clazz to select market data. If the domesticYCName is non-
        empty it is interpreted as indicating the name of the 'domstic'
        yield curve (see Model::getDomesticYCName()) */
    void getData(const IModel*                  model, 
                 const smartPtr<MarketData>&    market,
                 const CClassConstSP&           clazz,
                 const string&                  domesticYCName);

    /** Same as above, slightly different type for MarketData */
    void getData(const IModel*                  model, 
                 const MarketData*              market,
                 const CClassConstSP&           clazz,
                 const string&                  domesticYCName);

    /** Same as above, but no change in domesticYCName */
    void getData(const IModel*                  model, 
                 const smartPtr<MarketData>&    market,
                 const CClassConstSP&           clazz);

    /** Same as above, but no change in domesticYCName */
    void getData(const IModel*                  model, 
                 const MarketData*              market,
                 const CClassConstSP&           clazz);

    /** Returns true if the wrapper is getting its market data from the
        cache */
    bool usingCache() const;

    /** sets whether the wrapper should use the cache to retrieve the
        corresponding market data or not */
    void setCacheUse(bool useCache);

    /** Returns the name of the market object contained in this
        wrapper */
    string getName() const;

    /** Used for building instances of market object wrappers from
        strings. requiredType must be a MarketObject (or derived from
        a MarketObject) */
    static IObject* createFromString(CClassConstSP requiredType,
                                     const string& data);

    // is this an empty wrapper i.e. no name and no object
    bool isEmpty() const;

    //// Routes through equalTo. Method added to support
    //// instantiating array template
    bool operator==(const MarketObjectWrapper& rhs) const;

    /** We need it to sort the array of MarketObjects, 
        e.g. MarketWrapper<X>s */
    bool operator<(const MarketObjectWrapper & g2Compare) const;

protected:
    MarketObjectWrapper(const CClassConstSP& clazz);
    MarketObjectWrapper(const CClassConstSP& clazz, const string& name);
    MarketObjectWrapper(const CClassConstSP&  clazz, 
                        bool                  useCache);
    /** Throws an exception with message getClass()->getName()+ " with name "+
        getName()+" is null")*/
    void throwExceptionForNullObject() const;
private:
    MarketObjectWrapper();
    static void load(CClassSP& clazz);
    string          name;
    bool            useCache; // $unregistered
};

// typedef for smart pointers to MarketObjectWrapper
typedef smartConstPtr<MarketObjectWrapper> MarketObjectWrapperConstSP;
typedef smartPtr<MarketObjectWrapper> MarketObjectWrapperSP;

DRLIB_END_NAMESPACE

//// backward compatibility (contents previously lived here)
#include "edginc/MarketWrapper.hpp"

#endif
