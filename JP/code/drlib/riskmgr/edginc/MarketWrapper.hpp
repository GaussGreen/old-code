//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MarketWrapper.hpp
//
//   Description : template class for wrappers for derived instances of 
//                 MarketObjects
//
//   Author      : Mark A Robson
//
//   Date        : 19 Mar 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_MARKET_WRAPPER_HPP
#define EDG_MARKET_WRAPPER_HPP

#include "edginc/MarketObject.hpp"
#include "edginc/AtomicArray.hpp"

DRLIB_BEGIN_NAMESPACE
// template class for wrappers for derived instances of MarketObjects
template <class X> class MarketWrapper: public MarketObjectWrapper{
public:
    static CClassConstSP const TYPE;

    virtual ~MarketWrapper();

    /** Creates a default empty market object wrapper  */
    explicit MarketWrapper();

    /** Creates a market object wrapper where the market data is to be
        retrieved from the market data cache */
    explicit MarketWrapper(const string& name);

    /** Creates a market object wrapper where the market data is already
        within the wrapper - does not clone supplied object */
    explicit MarketWrapper(const smartPtr<X>& object);

    /** Creates a market object wrapper where the market data is already
        within the wrapper - takes ownership of the supplied object */
    explicit MarketWrapper(X* object);
    
    /** Creates a market object wrapper from existing wrapper
        - does not clone supplied object */
    MarketWrapper(const MarketWrapper& rhs);

    /** Returns the smart pointer within the wrapper. Does NOT throw
        an exception if null */
    virtual MarketObjectSP getMO() const;

    /** Stores the supplied object in the wrapper - no copy should be made */
    virtual void setObject(const MarketObjectSP& marketObject);

    /** Returns the smart pointer within the wrapper. Does NOT throw
        an exception if null */
    smartConstPtr<X> getSP() const;

    /** Returns the smart pointer within the wrapper. Does NOT throw
        an exception if null */
    smartPtr<X> getSP();

    /** Returns the MarketObject* within the wrapper. Throws an
        exception if MarketObject is NULL */
    const X* operator->() const;

    /** Returns the MarketObject* within the wrapper. Throws an
        exception if MarketObject is NULL */
    X* operator->();

    /** Returns the MarketObject* within the wrapper */
    const X* get() const;

    /** Returns the MarketObject* within the wrapper */
    X* get();

    /** Returns null if the market object inside is null */
    bool operator!() const;

    /** Due to strange C++ scoping rules, we need this statement to stop this
        method hiding the parent's methods */
    using MarketObjectWrapper::getData;

    /** Retrieves market data from cache if needed. Asks for data of type of
        supplied clazz to select market data */
    void getData(const IModel*              model, 
                 const MarketData*          market);

    /** Retrieves market data from cache if needed. Asks for data of type of
        supplied clazz to select market data */
    void getData(const IModel*                  model, 
                 const smartPtr<MarketData>&    market);

    /** Retrieves market data from cache if needed. Asks for data of type of
        supplied clazz to select market data. See MarketObject::getData
        for details on domesticYCName parameter */
    void getData(const IModel*              model, 
                 const MarketData*          market,
                 const string&              domesticYCName);

    /** Retrieves market data from cache if needed. Asks for data of type of
        supplied clazz to select market data */
    void getData(const IModel*                  model, 
                 const smartPtr<MarketData>&    market,
                 const string&                  domesticYCName);

    /** Returns the type of the market object that this is a wrapper for */
    virtual CClassConstSP getMOType() const;

private:
    smartPtr<X>  object;

    static IObject* defaultMarketWrapper();

    DLL_FIX_FOR_TEMPLATE_TYPE; // work around for VC71 and dlls
    static void load(CClassSP& clazz);

};
/** 
  * Returns all names of the market wrappers (same order) *********************
**/
template <class X>  
        StringArraySP getNames( const array< MarketWrapper < X> > & gWrappers);
/** 
  * Returns index of a market wrapper with name gName in the array ************
**/
template <class X> int getIndexOfNameInWrapperArray(
                                const string & gName, 
                                const array< MarketWrapper < X> > & gWrappers);


// then include the actual template implementation (allows us to use extern
// on template)
#include "edginc/MarketWrapper.inl"

DRLIB_END_NAMESPACE
#endif
