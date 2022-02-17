//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : BarrierUtil.hpp
//
//   Description : 
//
//   Date        : Q Hou Feb 2005
//
//
//----------------------------------------------------------------------------

#ifndef EDR_BarrierUtil_HPP
#define EDR_BarrierUtil_HPP

#include "edginc/Object.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/Array.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/EventResults.hpp"

DRLIB_BEGIN_NAMESPACE

/****************************************************************************/
// interface classes

class CashFlow;
class PhysicalDelivery;
class MCPathGenerator;

class PRODUCTS_DLL IBarrierUtil {
public:
    virtual ~IBarrierUtil(){};

    // relevant dates from this barrier itself
    virtual void relevantDates(DateTimeArray &relevantDts, bool &hasContinuous) const = 0;

    // linkage, memory management w.r.t. simulation timeline
    virtual void preprocess(DateTime valueDate, DateTimeArrayConstSP simDates) = 0;

    // return sim dates that's relevant to this barrier
    virtual const DateTimeArray& relevantMonDates() const = 0;

    // update path, may compute all required hitValues for a new path
    virtual void pathUpdated(const MCPathGenerator* pathGenIn) = 0;
    
    // get and/or calc hitValue
    virtual double hitValue(int obsStep) const = 0;

    virtual bool getHitDataForBreachEvent(double& refLevel, double& barrierLevel, double& assetlevel, string& monType) const = 0;

    virtual bool getEvents( EventResults* events, bool useIsHitFlags, const string& barrierName,
        StringArray* assNames, const DateTime& hitDate, int maxHits) const = 0;

    virtual bool getHitDateAndStep(int& step, DateTime& hitDate) const=0;

    class PRODUCTS_DLL ICanHaveBend {
    public:
        virtual ~ICanHaveBend(){};
        virtual void bendDates(DateTimeArray &startDates, DateTimeArray &endDates) const = 0;
        virtual double hitValueBend(int obsStep, int payObsStep) const = 0;
    };

    class PRODUCTS_DLL ICanHaveSticky {
    public:
        virtual ~ICanHaveSticky(){};
        virtual bool isSticky() const = 0;
        virtual double initialHitVal() const = 0;
        virtual int periodEndStep(int obsStep) const = 0;
    };

};

typedef refCountPtr<IBarrierUtil> IBarrierUtilSP;

class PRODUCTS_DLL IBarrierPay {
public:
    virtual ~IBarrierPay(){};
    
    virtual void preprocess(DateTime valueDate, 
                            DateTimeArrayConstSP simDates,
                            const DateTimeArray& monitorDates,
                            const YieldCurve* discount) = 0;
    
    // calc (discounted) value of payment associated with a step
    // return this value if payment in future
    // if payment/physDelivery is fixed due to historical events, return that
    virtual double value(int payObsStep,
        DateTime valueDate,
        const MCPathGenerator* pathGen) = 0;
    
    // collect knownCF and physical delivery
    virtual void collect(int obsStep, 
        const MCPathGenerator* pathGenIn,
        double scalingFactor,
        array<CashFlow, CashFlow>& knownCFs, 
        array<PhysicalDelivery>& phyDs) = 0;
    
    
    // query class to help barrier pay creator to create appropriate barrier and pay object
    virtual bool payOnce() const = 0;
    
    virtual bool payIfHitOne() const = 0; // true if payment occur when hitValue = 1
    
    // return dates when spot is needed to determine payments, not the
    // payment date or physical delivery date, which normally are later due to settlement
    virtual void payObsDates(DateTimeArray &obsDts, bool &hasContinuous) const = 0;

    // return monitor dates that's relevant for payment
    virtual const DateTimeArray& relevantMonDates() const = 0;
};

typedef refCountPtr<IBarrierPay> IBarrierPaySP;

class PRODUCTS_DLL IBarrierMgr {
public:
    virtual ~IBarrierMgr(){};

    virtual void preprocess(DateTime valueDate, DateTimeArrayConstSP simDates, const YieldCurve* discount) = 0;

    virtual double calculate(DateTime valueDate, const MCPathGenerator* pathGenIn) = 0;

    virtual void getKnowCFPhyD(DateTime valueDate, const MCPathGenerator* pathGenIn, double scalingFactor,
                               array<CashFlow, CashFlow>& knownCFs, array<PhysicalDelivery>& phyDs) = 0;

    virtual void getEvents(const  MCPathGenerator*  pathGen,
        DateTime valueDate,
        EventResults* events,
        const string& barrierName,
        StringArray* assNames) = 0;
};

typedef refCountPtr<IBarrierMgr> IBarrierMgrSP;


class PRODUCTS_DLL BarrierFactory {
public:
    // make some building blocks
    static IBarrierUtilSP makeSticky(IBarrierUtilSP barrier);
    static IBarrierUtilSP makeCoupler(vector<IBarrierUtilSP> barriers, int numHit, bool isOut, bool isUp);
    static IBarrierUtilSP makeReverser(IBarrierUtilSP barrier);

    // trivial pay give hitValue=1 for the requested dates
    // hence the manager will count/sum the hit values of the dates
    // which are the last prev mon date for each obs date
    static IBarrierPaySP makeTrivialPay(const DateTimeArray& obsDts);

    static IBarrierMgrSP makeManager(IBarrierUtilSP barrier, IBarrierPaySP pay);
};

DRLIB_END_NAMESPACE
                     
#endif
                     
