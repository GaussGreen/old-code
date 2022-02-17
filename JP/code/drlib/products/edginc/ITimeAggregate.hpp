//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ITimeAggregate.hpp
//
//   Description : Captures collapse of TIME dimension from array to single value in various ways
//
//   Date        : Mar 04
//
//
//----------------------------------------------------------------------------

#ifndef EDR_ITIMEAGGREGATE_HPP
#define EDR_ITIMEAGGREGATE_HPP

#include "edginc/Object.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/InstrumentSettlement.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/IDoubleArray.hpp"

DRLIB_BEGIN_NAMESPACE

// View of IMCProduct which provides required info
class PRODUCTS_DLL ITimeAggProductView {
public:    
    /** Returns the value date ie today */
    virtual const DateTime& getValueDate() const = 0;
    
    // Could offer a pv() method only... XXX
    virtual const YieldCurve* getYieldCurve() const = 0;
    
    virtual DateTime settles(const DateTime& aDate) const = 0;

};

// Would like this to be derived from IAggregate, but maybe that wouldn't work...
//class ITimeAggregate : virtual public IAggregate {
class PRODUCTS_DLL ITimeAggregate {
public:
    virtual double aggregate() = 0;

    // To support the PAYMENT_DATES request
    virtual const DateTimeArray* getPaymentDates() const = 0;

    // To support the KNOWN_CASHFLOWS request
    virtual const CashFlowArray* getKnownFlows() const = 0;

    virtual ~ITimeAggregate() {};
};

typedef refCountPtr<ITimeAggregate> ITimeAggregateSP;

// Build ITimeAggregate via a Maker class so the ITimeAggregate can
// be identified once and aggregate() can act on known memory. 
// We're building the thing we need in 2 stages.
// class ITimeAggregateMaker : virtual public IAggregateMaker {
class PRODUCTS_DLL ITimeAggregateMaker : virtual public IObject {
public:
    static CClassConstSP const TYPE;
    
    // The simulation will need to know when to place dates
    virtual const DateTimeArray& getDates() const = 0;

    virtual ITimeAggregate* getAggregate(const DateTime&            valueAtDate,
                                         const ITimeAggProductView* prodView,
                                         IDoubleArray*              comps) = 0;

    virtual ~ITimeAggregateMaker() {};

private:
    static void load(CClassSP& clazz);
};

typedef smartPtr<ITimeAggregateMaker> ITimeAggregateMakerSP;

/*********************************************************************/
// For values corresponding to each date, sum the PVs
// Allow "payAtMat" choice.
// Note that in order to allow subsequent processing to interpret values correctly
// I'm going to have the aggregate valued at the final date (e.g. any overall option
// will apply a strike at that date, and to have this value today will mean effectively
// PVing the strike)
// XXX Perhaps should have a "valueAtDate" for all TimeAggs (passed in getAggregate?) which 
// XXX defines where ta->aggregate()'s value applies.
class PRODUCTS_DLL SumPVAggregateMaker : public CObject,
                            virtual public ITimeAggregateMaker {
public:
    static CClassConstSP const TYPE;
    friend class SumPVAggregate;
    friend class SumPVAggregateMakerHelper; // for hiding the registration

    void validatePop2Object();

    virtual const DateTimeArray& getDates() const;

    virtual ITimeAggregate* getAggregate(const DateTime&            valueAtDate,
                                         const ITimeAggProductView* prodView,
                                         IDoubleArray*              comps);

protected:
    DateTimeArray aggDates;
    DateTimeArray aggPayDates;
    bool          adjustPayDatesForInstSettlement;

private:
    SumPVAggregateMaker(): CObject(TYPE), aggDates(0), aggPayDates(0),
        adjustPayDatesForInstSettlement(true) {} // for reflection
    SumPVAggregateMaker(const SumPVAggregateMaker& rhs); // not implemented
    SumPVAggregateMaker& operator=(const SumPVAggregateMaker& rhs); // not implemented

};

typedef smartPtr<SumPVAggregateMaker> SumPVAggregateMakerSP;
typedef array<SumPVAggregateMakerSP, SumPVAggregateMaker> SumPVAggregateMakerArray;
typedef smartPtr<SumPVAggregateMakerArray> SumPVAggregateMakerArraySP;

DRLIB_END_NAMESPACE

#endif





