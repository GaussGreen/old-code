//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MarketDataFetcher.hpp
//
//   Description : Helper class for models to get data out of market cache
//
//   Author      : Andrew J Swain
//
//   Date        : 1 February 2002
//
//----------------------------------------------------------------------------

#ifndef MARKETDATAFETCHER_HPP
#define MARKETDATAFETCHER_HPP

#include "edginc/VirtualDestructorBase.hpp"
#include "edginc/MarketObject.hpp"
#include "edginc/CorrelationBase.hpp"
#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE

class MarketData;
class IModel;

/**
Introduction:
The primary purpose of MarketDataFetchers is encapsulated by the 'fetch' method.
In particular, clients ask for a piece of market data of a certain type to be
retrieved. The MarketDataFetchers job is to decide whether or not to fetch the
market data (if not, null is returned) and if so then whether a more specialised
type should be requested. For example a CAsset typically wants a CVolBase and
a MarketDataFetcher might need to select a particular [derived] type of CVolBase
since there may be multiple instances of objects in the market data cache which
have the same name and are derived from CVolBase.

Note that there is a MDFUtul class (in market) which has some utility methods
for common configurations.

Motivation:
The original implementation started off simple but increasingly got more 
complex and non-scalable. The approach here uses a 'declarative' approach
and clients can say get this type or that type of market data as well as 
specifying what types to specialise to.

Detail:
The internals of the MarketDataFetcher class has been refactored to allow
its use to be controlled in a declarative fashion. This means that at
construction time a series of rules are issued to the MarketDataFetcher. These
rules controls its behaviour when it is used to retrieve market data.
It is still possible to override the implementation in derived classes but
this should be needed infrequently.

Although there are two flavours of rules, they all essentially say the same
thing which is when given a type of market data (eg CVol) whether to get that
piece of market data or not, and if the market data is to be retrieved then what
type of market data to retrieve.
The two flavours differ in that one flavour allows the rule only to be used
for a component of market data when inside another piece of market data of
a specified type. This refers to the common situation when market data refers
to another piece of market data. When 'fetching' this latter piece of market
data you are 'inside' the former piece of market data.

When a MarketDataFetcher is created it picks up a set of default rules. Thus there are two ways of applying the rules. The first is to change the default rules that all MarketDataFetcher pick up when they are created, whilst the second changes an existing instance of a MarketDataFetcher.

The rules for the first flavour are configured using the methods

        Select whether market data of the specified type should be, by
        default, retrieved or not from the MarketData cache and, if it is,
        what specialised type should be retrieved. Use NULL for typeToGet
        to indicate that the objects of requestedType should have the type
        requested unchanged. This is a global setting and this method can
        only be called at startup (ie during load methods).  It
        [potentially] affects all market data fetchers so care needs to be
        exercised (otherwise you may break every test case ...)
    static void setDefaultRetrievalMode(CClassConstSP requestedType,
                                        bool          retrieve,
                                        CClassConstSP typeToGet);

        Select whether market data of the specified type should be
        retrieved or not from the MarketData cache by this
        MarketDataFetcher. This is the same as setDefaultRetrievalMode but is
        specific to this instance of the MarketDataFetcher
    void setRetrievalMode(CClassConstSP requestedType, 
                          bool          retrieve,
                          CClassConstSP typeToGet) const; // hacky 

Whilst the rules for the second flavour are configured using the methods:

        This is the same as the above setDefaultRetrievalMode
        method but the rule for
        retrieving or not retrieving the object is only used if, when
        retrieving data, we are inside an object of the specified
        containing type. For example, this is used to avoid getting
        IRVols inside yield curves. This is a global setting and this
        method can only be called at startup (ie during load methods).
        It [potentially] affects all market data fetchers so care
        needs to be exercised (otherwise you may break every test case
        ...)
    static void setDefaultRetrievalMode(CClassConstSP requestedType,
                                        CClassConstSP containingType,
                                        bool          retrieve,
                                        CClassConstSP typeToGet);

        This is the same as the above getRetrievalMode method but the
        rule for retrieving or not retrieving the object is only used
        if, when retrieving data, we are inside an object of the
        specified containing type. For example, this is used to avoid
        getting IRVols inside yield curves
    void setRetrievalMode(CClassConstSP requestedType,
                          CClassConstSP containingType,
                          bool          retrieve,
                          CClassConstSP typeToGet) const; //


Subtleties
There are some subtleties in the way the rules interact with each
other. The order below also reflects the order in which the relevant rule is
determined.

1. It is possible to have multiple matching rules for a given type. In such a situation the rule for the 'closest' type is used. Here if type A derives from P1 and P2 then P1 is closest to A if P1 derives from P2, whilst P2 is closest to A if P2 derives from P1. If P1 is not derived from P2 or
P1 is not derived from P2 then an error is thrown.

2. When selecting whether or not to retrieve an item of market data the market
data fetcher may be configured to choose to choose a particular type. For
example, it may be configured to get a VolHyperTrig vol when a CVolBase (or a type derived from CVolBase) is asked for. However, this rule will be ignored if the actual incoming type that is requested is not compatible with the type it is
to be switched too. (Continuing the example, if a request for a VolSurface is made then the rule to switch CVolBase to VolHyperTrig will be ignored since a VolHyperTrig cannot be used where a VolSurface is required.)

3. The rules which take affect only 'inside' certain objects have precedence over the general rules. For example if a MarketDataFetcher is configured to get IR swaption vols in general but also not to get IRVolBases inside YieldCurves then the latter rules takes precedence inside yield curves.

4. If there are multiple matching rules of the type described in (3) then the
rule for the outermost object takes precedence (this seemed more useful than the other way round and is used by the MarketDataFetcherCDS for quanto adjusted curves). 

5. If no rule is found, then an attempt is made to get the requested piece of
market data of the original requested type.

Examples
--------
1. To get/not get CurrencyBasis objects
setRetrievalMode(CurrencyBasis::TYPE, getCurrencyBasis, NULL);
2. To set whether or not IRVols within yield curves are retrieved
    setRetrievalMode(IRVolBase::TYPE,
                     IYieldCurve::TYPE,
                     getStochasticYieldCurves, NULL);
3. See MarketDataFetcherSRM::setSwaptionVolFlag(bool getSwaptionVols)
4. See MDFUtil.cpp (in market)
 */
class RISKMGR_DLL MarketDataFetcher: public virtual VirtualDestructorBase {

public:
    virtual ~MarketDataFetcher();

    /** Essentially a wrapper around MarketData::GetData(model, name, type).
        Retrieves the market object with
        the specififed name and type (and retrieves its market data). The
        domesticYCName flags what the current 'domestic' currency is (see
        getDomesticYCName() for more information). */
    virtual MarketObjectSP getMarketObject(
        const IModel*     model,
        const MarketData* market,
        const string&     name,
        CClassConstSP     type,
        const string&     domesticYCName) const;

    /** default fetch method. Note that implementation now uses a declarative
        approach to determine if and what to retrieve. The MarketDataFetcher
        can be configured using the setRetrievalMode methods below. */
    virtual MarketObjectSP fetch(const MarketData*    market,
                                 const string&        name,
                                 const CClassConstSP& type,
                                 const IModel*        model) const;

    /** By default, the same method on Model delegates to this method.
        Briefly, it allows the MarketDataFetcher context information 
        regarding what object it is currently in */
    virtual void getComponentMarketData(const IModel*        model,
                                        const MarketData*    market,
                                        MarketObjectSP       mo) const;

    /** Same as above but specifies a 'domestic' yield curve for the 
        current component (see getDomesticYCName()) */
    virtual void getComponentMarketData(
        const IModel*     model,
        const MarketData* market,
        MarketObjectSP    mo,
        const string&     domesticYCName) const;

    /** Select whether market data of the specified type should be, by
        default, retrieved or not from the MarketData cache and, if it is,
        what specialised type should be retrieved. Use NULL for typeToGet
        to indicate that the objects of requestedType should have the type
        requested unchanged. This is a global setting and this method can
        only be called at startup (ie during load methods).  It
        [potentially] affects all market data fetchers so care needs to be
        exercised (otherwise you may break every test case ...) */
    static void setDefaultRetrievalMode(CClassConstSP requestedType,
                                        bool          retrieve,
                                        CClassConstSP typeToGet);

    /** This is the same as the above method but the rule for
        retrieving or not retrieving the object is only used if, when
        retrieving data, we are inside an object of the specified
        containing type. For example, this is used to avoid getting
        IRVols inside yield curves. This is a global setting and this
        method can only be called at startup (ie during load methods).
        It [potentially] affects all market data fetchers so care
        needs to be exercised (otherwise you may break every test case
        ...)  */
    static void setDefaultRetrievalMode(CClassConstSP requestedType,
                                        CClassConstSP containingType,
                                        bool          retrieve,
                                        CClassConstSP typeToGet);

    /** Select whether market data of the specified type should be
        retrieved or not from the MarketData cache by this
        MarketDataFetcher. This is the same as setDefaultRetrievalMode but is
        specific to this instance of the MarketDataFetcher  */
    void setRetrievalMode(CClassConstSP requestedType, 
                          bool          retrieve,
                          CClassConstSP typeToGet) const /* hacky */;

    /** Query whether market data of the specified type is retrieved or
        not from the MarketData cache by this MarketDataFetcher ignoring
        any settings dependent on what object we are in when retrieving
        market data. This is the partner to setRetrievalMode method which
        takes 3 arguments. */
    bool getRetrievalMode(CClassConstSP  requestedType) const;

    /** This is the same as the above getRetrievalMode method but the
        rule for retrieving or not retrieving the object is only used
        if, when retrieving data, we are inside an object of the
        specified containing type. For example, this is used to avoid
        getting IRVols inside yield curves  */
    void setRetrievalMode(CClassConstSP requestedType,
                          CClassConstSP containingType,
                          bool          retrieve,
                          CClassConstSP typeToGet) const /* hacky */;

    /** Query whether market data of the specified type is retrieved or
        not from the MarketData cache by this MarketDataFetcher when
        inside an object of type containingType. This is the partner to
        setRetrievalMode method which takes 4 arguments. */
    bool getRetrievalMode(CClassConstSP  requestedType,
                          CClassConstSP  containingType) const;
    /** Sets whether correl swap basis shall be used or not */
    virtual DateTimeSP setCorrSwapExpiry(DateTimeSP corrSwapExpiryInput) const;
    virtual DateTimeSP getCorrSwapExpiry() const;

    /** Return the name of the 'domestic' yield curve, if available, or
     *  an empty string if not available. */
    virtual const string& getDomesticYCName() const;

    /** Sets the 'domestic' yield curve name */
    virtual void setDomesticYCName(string domYCName) const;

    /** Retrieve a correlation from the cache with the supplied name and
     *  type specified by corrType. type1 and type2 are the types of the
     *  two objects with respect to which this is a correlation. The
     *  return correlation will be correcly configured for the
     *  appropriate sensitivity (eg PHI, FX_PHI) etc. Note that a model
     *  may decide that it does not use this type of correlation (eg IR v EQ) 
     *  and return a 'dummy' version if say the object does not live in the
     *  cache. In a similar manner, the returned object may be sensitive to
     *  no sensitivities if the model is not going to use the correlation */
    virtual CorrelationBaseSP getCorrelation(const IModel*     model,
                                             const string&     corrName,
                                             CClassConstSP     type1,
                                             CClassConstSP     type2,
                                             CClassConstSP     corrType,
                                             const MarketData* market) const;

    /** Invoked for each piece of market data (whether already inline in
     *  instrument or pulled from cache). The default implementation just
     *  returns mo. Derived classes can use this to replace market data
     *  objects with other instances. It compliments GetMarket in that this
     *  method works for instruments which have the market data inside them
     *  already */
    virtual MarketObjectSP modifyMarketData(
                                         const IModel*         model,
                                         const MarketData*     market,
                                         const CClassConstSP&  clazz,    
                                         const MarketObjectSP& mo) const;
    //// uses configuration as specified by setDefaultRetrievalMode invocations
    MarketDataFetcher(); // defaults allowVarSwapBasis and getStochasticYCs to false

    /** Resolve an array of multiply selected market data
        into a single chosen market datum */
    virtual MarketObjectConstSP resolveMultiData(const MarketObjectArray& marketObjs,
                                                 CClassConstSP            type) const;

protected:
    /** Given the requested type, works out what type to actually get. A null
        return value means not to retrieve the item at all */
    CClassConstSP getTypeToRetrieve(CClassConstSP requestedType) const;

    //// this vanishes once get rid of StochasticYieldCurves, but until it does
    //// you need to call this before any calls to getTypeToRetrieve 
    CClassConstSP stochasticYCFix(
        const MarketData*    market,
        const string&        name,
        CClassConstSP        type) const;
private:
    MarketDataFetcher(const MarketDataFetcher& rhs);
    MarketDataFetcher& operator=(const MarketDataFetcher& rhs);
    class RetrievalConfig;
    class Imp;
    friend class Imp;
    auto_ptr<Imp> my; // hides implementation
};

typedef smartConstPtr<MarketDataFetcher> MarketDataFetcherConstSP;
typedef smartPtr<MarketDataFetcher> MarketDataFetcherSP;
#ifndef QLIB_MARKETDATAFETCHER_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartPtr<MarketDataFetcher>);
#else
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartPtr<MarketDataFetcher>);
#endif
DRLIB_END_NAMESPACE
#endif
