//----------------------------------------------------------------------------
//
//   File        : IForwardRatePricer.hpp
//
//   Description : A general interface for models that are capable
//                 of generating and pricing forward rates
//                 Also includes IHasForwardRatePricer interface
//                 that pricing models can implement if they contain
//                 an IForwardRatePricer, which makes it accessible
//                 Also includes IWantsForwardRatePricer interface
//                 that allows a forward rate object to declare the
//                 type of model it requires
//
//----------------------------------------------------------------------------

#ifndef QLIB_IFORWARDRATEPRICER_HPP
#define QLIB_IFORWARDRATEPRICER_HPP

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/DECLARE.hpp"
#include "edginc/Object.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(MarketData);
FORWARD_DECLARE(IModel);

//-------------------
// IForwardRatePricer
//-------------------

class RISKMGR_DLL IForwardRatePricer : virtual public IObject {
public:
    static CClassConstSP const TYPE;

    static const string INTEREST_RATE_FORWARD;
    static const string CREDIT_FORWARD;
    static const string CLOSED_FORM_FORWARD;

    IForwardRatePricer();
    virtual ~IForwardRatePricer();

    virtual void Forward(CObject *instrument, double *result) = 0;
    virtual void getMarket(IModel *model, const MarketData *market);


private:
    static void load(CClassSP& clazz);
};

DECLARE(IForwardRatePricer)

//----------------------
// IHasForwardRatePricer
//----------------------

class RISKMGR_DLL IHasForwardRatePricer : virtual public IObject {
public:
    static CClassConstSP const TYPE;

    IHasForwardRatePricer();
    virtual ~IHasForwardRatePricer();

    /** Key method providing access to the pricer */
    virtual IForwardRatePricerSP getForwardRatePricer() const = 0;

protected:
    /** implementing classes can make use of this
     *  utility method to instantiate a default pricer
     *  Not expected to be overridden
     */
    IForwardRatePricerSP getDefaultForwardRatePricer() const;

private:
    static void load(CClassSP& clazz);
};

//------------------------
// IWantsForwardRatePricer
//------------------------

class RISKMGR_DLL IWantsForwardRatePricer : virtual public IObject {
public:
    static CClassConstSP const TYPE;

    IWantsForwardRatePricer();
    virtual ~IWantsForwardRatePricer();

    /** Key method allowing the object to specify the
      * model it requires */
    virtual string requiredPricer() const = 0;

private:
    static void load(CClassSP& clazz);
};

DRLIB_END_NAMESPACE

#endif
