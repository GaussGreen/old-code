//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IRebate.hpp
//
//   Description : 
//
//   Date        : Oct 2003
//
//
//----------------------------------------------------------------------------

#ifndef EDR_IREBATE_HPP
#define EDR_IREBATE_HPP

#include "edginc/Object.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/Schedule.hpp"
#include "edginc/SwapLegIntFace.hpp"

DRLIB_BEGIN_NAMESPACE

class IBarrierPay;

/*************************************************************/
// Allow some flexibility in styles of rebate by providing interface and 
// opportunity to extend choice of implementers. Need to consider
// hitting time, different per asset etc.
class PRODUCTS_DLL IRebate {
public:
    virtual double getLevel(int simDateIdx) const = 0;

    // if hit, when pay. NOT settle adjusted
    virtual DateTime getPayDate(const DateTime& hitDate) const = 0;

    virtual ~IRebate() {};
};

typedef refCountPtr<IRebate> IRebateSP;

// Build IRebate via a Maker class.
// It's a bit like building the thing we need in 2 stages.
class PRODUCTS_DLL IRebateMaker : virtual public IObject {
public:
    static CClassConstSP const TYPE;
    
    virtual IRebateSP getRebate(const DateTimeArray& simDates,
                                const YieldCurve*    discount) = 0;

    // populate from market cache 
    virtual void getMarket(const IModel*     model, 
                           const MarketData* market) = 0;

    // support barrier pay
    static refCountPtr<IBarrierPay> createBarrierPay(const IRebateMaker* rebateMaker, 
                                                     bool isOut, 
                                                     const InstrumentSettlement* instSettle);

    virtual ~IRebateMaker() {};

private:
    static void load(CClassSP& clazz);
};

typedef smartPtr<IRebateMaker> IRebateMakerSP;

class PRODUCTS_DLL FlatRebateMaker : public CObject,
                        virtual public IRebateMaker {
public:
    static CClassConstSP const TYPE;
    friend class FlatRebate;
    friend class FlatRebateMakerHelper;

    FlatRebateMaker(double rebateAmount);

    virtual IRebateSP getRebate(const DateTimeArray& simDates,
                                const YieldCurve*    discount);

    virtual void getMarket(const IModel*     model, 
                           const MarketData* market);

    double getRebateAmount() const;

protected:
    double rebateAmount;
    bool   rebatePayAtHit; // $unregistered

    FlatRebateMaker(CClassConstSP clazz);

private:
    FlatRebateMaker(): CObject(TYPE), rebateAmount(0.0), rebatePayAtHit(true) {} // for reflection
    FlatRebateMaker(const FlatRebateMaker& rhs); // not implemented
    FlatRebateMaker& operator=(const FlatRebateMaker& rhs); // not implemented
};

typedef smartPtr<FlatRebateMaker> FlatRebateMakerSP;

class PRODUCTS_DLL ScheduleRebateMaker : public CObject,
                            virtual public IRebateMaker {
public:
    static CClassConstSP const TYPE;
    friend class ScheduleRebate;
    friend class ScheduleRebateMakerHelper;

    virtual void validatePop2Object();

    virtual IRebateSP getRebate(const DateTimeArray& simDates,
                                const YieldCurve*    discount);

    virtual void getMarket(const IModel*     model, 
                           const MarketData* market);
    
    // make scheduleSP
    ScheduleSP makeSchedule() const;

protected:
    DateTimeArray rebateSchedDates;
    DoubleArray   rebateSchedValues;
    string        rebateSchedInterp;

    // optional field for discrete rebate pay period dates
    // hits within a period are paid at end of period
    // empty date indicate payAtHit!
    DateTimeArray rebatePayPeriodDates;

    ScheduleRebateMaker(CClassConstSP clazz);

private:
    ScheduleRebateMaker(): CObject(TYPE) {} // for reflection
    ScheduleRebateMaker(const ScheduleRebateMaker& rhs); // not implemented
    ScheduleRebateMaker& operator=(const ScheduleRebateMaker& rhs); // not implemented

};

typedef smartPtr<ScheduleRebateMaker> ScheduleRebateMakerSP;


/*********************************************************************/
/** Rebate is "FloatStream + FixedStream + Fixed Amount"             */
/*********************************************************************/
class PRODUCTS_DLL StreamRebateMaker : public CObject,
                          virtual public IRebateMaker {
public:
    static CClassConstSP const TYPE;
    friend class StreamRebate;
    friend class StreamRebateMakerHelper;
        
    virtual IRebateSP getRebate(const DateTimeArray& simDates,
                                const YieldCurve*    discount);

    virtual void getMarket(const IModel*     model, 
                           const MarketData* market);

protected:
    // floating leg
    LiborLegSP        floatStreamRebate;
    string            floatKoStubRule;
    // fixed leg
    FixedLegSP        fixedStreamRebate;
    string            fixedKoStubRule;
    // fixed amount
    DateTimeArray rebateSchedDates;
    DoubleArray   rebateSchedValues;
    string        rebateSchedInterp;

private:
    StreamRebateMaker(); // for reflection
    StreamRebateMaker(const StreamRebateMaker& rhs); // not implemented
    StreamRebateMaker& operator=(const StreamRebateMaker& rhs); // not implemented
};

typedef smartPtr<StreamRebateMaker> StreamRebateMakerSP;

DRLIB_END_NAMESPACE

#endif

