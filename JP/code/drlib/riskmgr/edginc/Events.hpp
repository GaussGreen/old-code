
#ifndef EDG_EVENTS_H
#define EDG_EVENTS_H

#include "edginc/Object.hpp"
#include "edginc/IDoubleArray.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/CashFlow.hpp"

DRLIB_BEGIN_NAMESPACE

class IModel;
class CInstrument;
class EventResults;

// interface actual events must satisfy 
class RISKMGR_DLL IEvent: public virtual IObject {
public:
    static CClassConstSP const TYPE;

    // must implement this to define specific event interface
    virtual CClassConstSP getEventInterface() const = 0;

    // gets the event for this interface
    virtual void getEvents(IObjectConstSP       obj,
                           IModel*              model,
                           const DateTime&      eDate,
                           EventResults*        events) const = 0;

    // if we get risk and legal events we need to know whether a risk
    // event will just be a duplicated legal event (and hence removable)
    virtual bool isValidRiskEvent(const DateTime& today) const = 0;
};

/*****************************************************************************/
// Base class for various flavours of MATRIX event
class RISKMGR_DLL Event: public CObject{
public:
    static CClassConstSP const TYPE;

protected:
    // for reflection
    Event(CClassConstSP clazz);

    //for inheritance
    Event(CClassConstSP clazz, const DateTime& eDate);
    DateTime    eventDate;

private:
    static void load(CClassSP& clazz);
    Event(const Event &rhs);
    Event& operator=(const Event& rhs);

    static IObject* defaultEvent();
};

typedef smartPtr<Event> EventSP;
typedef array<EventSP, Event> EventArray;

/******************************************************************************/
class RISKMGR_DLL BarrierBreach: public Event, public virtual IEvent{
public:
    static CClassConstSP const TYPE;
    static const string KNOCK_IN;
    static const string KNOCK_OUT;
    static const string LOCK_IN;
    static const string NOT_APPLICABLE;
    static const string DAILY;
    static const string EUROPEAN;
    static const string CONTINUOUS;

    static void load(CClassSP& clazz);

    // interface instrument must satisfy to handle this event
    class RISKMGR_DLL IEventHandler {
    public:
        static CClassConstSP const TYPE;

        virtual void getEvents(const BarrierBreach*     breach,
                               IModel*                  model,
                               const DateTime&          eDate,
                               EventResults*            events) const = 0;
    };

    BarrierBreach(const DateTime& eDate, string barrDesc, string monType,
                  string barrType, bool isItUp, int hitsLeft, StringArraySP names,
                  DoubleArraySP levels, DoubleArraySP bLevels);

    CClassConstSP getEventInterface() const;

    virtual void getEvents(IObjectConstSP       obj,
                           IModel*              model,
                           const DateTime&      eDate,
                           EventResults*        events) const;

    // if we get risk and legal events we need to know whether a risk
    // event will just be a duplicated legal event (and hence removable)
    virtual bool isValidRiskEvent(const DateTime& today) const;

private:
    BarrierBreach();
    BarrierBreach(const BarrierBreach &rhs);
    BarrierBreach& operator=(const BarrierBreach& rhs);
    static IObject* defaultBarrierBreach();

    // general data about the barrier 'breached'
    string          barrDesc;
    string          barrierType;
    string          monType;
    bool            isUp;
    int             hitsRemaining;

    // data about the assets that breached
    StringArraySP   names;
    DoubleArraySP   levels;
    DoubleArraySP   barrierLevels;
};

/******************************************************************************/
class RISKMGR_DLL FlexBarrierBreach: public Event, public virtual IEvent{
public:
    static CClassConstSP const TYPE;

    FlexBarrierBreach(const DateTime& eDate, string varName, double barrLevel,
                        double varLevel, bool isItUp);

    // interface instrument must satisfy to handle this event
    class RISKMGR_DLL IEventHandler {
    public:
        static CClassConstSP const TYPE;

        virtual void getEvents(const FlexBarrierBreach* breach,
                               IModel*                  model,
                               const DateTime&          eDate,
                               EventResults*            events) const = 0;
    };

    CClassConstSP getEventInterface() const;

    virtual void getEvents(IObjectConstSP       obj,
                           IModel*              model,
                           const DateTime&      eDate,
                           EventResults*        events) const;

    // if we get risk and legal events we need to know whether a risk
    // event will just be a duplicated legal event (and hence removable)
    virtual bool isValidRiskEvent(const DateTime& today) const;

private:
    FlexBarrierBreach();
    FlexBarrierBreach(const FlexBarrierBreach &rhs);
    FlexBarrierBreach& operator=(const FlexBarrierBreach& rhs);
    static IObject* defaultFlexBarrierBreach();
    static void load(CClassSP& clazz);

    string          variableName;
    bool            isUp;
    double          variableLevel;
    double          barrierLevel;
};

/******************************************************************************/
class RISKMGR_DLL TargetRedemption: public Event, public virtual IEvent{
public:
    static CClassConstSP const TYPE;
    static const string KNOCK_IN;
    static const string KNOCK_OUT;
    static const string NOT_APPLICABLE;

    TargetRedemption(const DateTime& eDate, double cpn, double sumCpn,
                        double target, double bonusCpn, double redemptionVal = 0.0,
                        string floatingLegType = NOT_APPLICABLE);

    // interface instrument must satisfy to handle this event
    class RISKMGR_DLL IEventHandler {
    public:
        static CClassConstSP const TYPE;

        virtual void getEvents(const TargetRedemption*  tarn,
                               IModel*                  model,
                               const DateTime&          eDate,
                               EventResults*            events) const = 0;
    };

    CClassConstSP getEventInterface() const;

    virtual void getEvents(IObjectConstSP       obj,
                           IModel*              model,
                           const DateTime&      eDate,
                           EventResults*        events) const;

    // if we get risk and legal events we need to know whether a risk
    // event will just be a duplicated legal event (and hence removable)
    virtual bool isValidRiskEvent(const DateTime& today) const;

private:
    TargetRedemption();
    TargetRedemption(const TargetRedemption &rhs);
    TargetRedemption& operator=(const TargetRedemption& rhs);
    static IObject* defaultTargetRedemption();
    static void load(CClassSP& clazz);

    double          coupon;
    double          sumCoupons;
    double          targetLevel;
    double          bonus;
    double          redemption;
    string          floatingLeg;
};

/******************************************************************************/
class RISKMGR_DLL Callability : public Event, public virtual IEvent{
public:
    static CClassConstSP const TYPE;
    static const string CALLABLE;
    static const string PUTTABLE;
    static const string CALL;
    static const string PUT;

    Callability(const DateTime& notifDate, const DateTime& callDate, 
                string type, double amount);

    // interface instrument must satisfy to handle this event
    class RISKMGR_DLL IEventHandler {
    public:
        static CClassConstSP const TYPE;

        virtual void getEvents(const Callability*       call,
                               IModel*                  model,
                               const DateTime&          eDate,
                               EventResults*            events) const = 0;
    };

    CClassConstSP getEventInterface() const;

    virtual void getEvents(IObjectConstSP       obj,
                           IModel*              model,
                           const DateTime&      eDate,
                           EventResults*        events) const;

    // if we get risk and legal events we need to know whether a risk
    // event will just be a duplicated legal event (and hence removable)
    virtual bool isValidRiskEvent(const DateTime& today) const;

private:
    Callability();
    Callability(const Callability &rhs);
    Callability& operator=(const Callability& rhs);
    static IObject* defaultCallability();
    static void load(CClassSP& clazz);

    DateTime        notifDate;
    DateTime        callDate;
    string          type;
    double          amount;
};

/******************************************************************************/
class RISKMGR_DLL SPIFixing : public Event, public virtual IEvent{
public:
    static CClassConstSP const TYPE;
    static const string ZC_BOND_FIXING;
    static const string MM_RATE_FIXING;
    static const string BOND_FLOOR_FIXING;

    SPIFixing(const DateTime& eDate, string type, double amount);

    // interface instrument must satisfy to handle this event
    class RISKMGR_DLL IEventHandler {
    public:
        static CClassConstSP const TYPE;

        virtual void getEvents(const SPIFixing*         fixing,
                               IModel*                  model,
                               const DateTime&          eDate,
                               EventResults*            events) const = 0;
    };

    CClassConstSP getEventInterface() const;

    virtual void getEvents(IObjectConstSP       obj,
                           IModel*              model,
                           const DateTime&      eDate,
                           EventResults*        events) const;

    // if we get risk and legal events we need to know whether a risk
    // event will just be a duplicated legal event (and hence removable)
    virtual bool isValidRiskEvent(const DateTime& today) const;

private:
    SPIFixing();
    SPIFixing(const SPIFixing &rhs);
    SPIFixing& operator=(const SPIFixing& rhs);
    static IObject* defaultSPIFixing();
    static void load(CClassSP& clazz);

    string          type;
    double          amount;
};

/******************************************************************************/

class CashflowInfo;
typedef smartPtr<CashflowInfo> CashflowInfoSP;
typedef smartConstPtr<CashflowInfo> CashflowInfoConstSP;
typedef array<CashflowInfoSP, CashflowInfo> CashflowInfoArray;
typedef smartPtr<CashflowInfoArray> CashflowInfoArraySP;

class RISKMGR_DLL CashflowInfo : public CObject {
public:
    static CClassConstSP const TYPE;

    struct CfType { 
        enum Enum { UNSET, PRINCIPAL, COUPON, REDEEMER };
    };
    typedef BoxedEnum<CfType::Enum> CfTypeBoxedEnum;

    struct AmountType {
        enum Enum { UNSET, KNOWN, ESTIMATED, UNKNOWN };
    };
    typedef BoxedEnum<AmountType::Enum> AmountTypeBoxedEnum;

    string componentName;
    CfType::Enum cfType;
    double amount;
    DateTime date;
    AmountType::Enum amountType;
    double multiplierUsed;

    StringArray keys;
	StringArray values;

    void setMultiplier(double multiplier);
    CashFlow getCashFlow() {return CashFlow(date, amount);}

    // Downgrades amountType to newAmountType
    void updateAmountType(AmountType::Enum newAmountType);

    // Merges the (key/values, amountType) of "cfi" with the current instance
    // if cfi.amountType != UNKNOWN, otherwise just set amountType to UNKNOWN
    void merge(const CashflowInfo &cfi, string const &prefix); 

	void push(const string &key, const string &value);
	void push(const string &key, double value);
	void push(const string &key, bool value);
	void push(const string &key, DateTime value);
	CashflowInfo(const CClassConstSP &type=TYPE);

private:
	static IObject* defaultCashflowInfo(void);
	static void load(CClassSP& clazz);
};

/*************************/

class RISKMGR_DLL KnownCashflows : public Event, public virtual IEvent{
public:
    static CClassConstSP const TYPE;

    KnownCashflows(const DateTime& eDate, CashFlowArraySP cashFlows, 
        string ccy, CashflowInfoArraySP info=CashflowInfoArraySP(   ));

    // interface instrument must satisfy to handle this event
    class RISKMGR_DLL IEventHandler {
    public:
        static CClassConstSP const TYPE;

        virtual void getEvents(const KnownCashflows*    flows,
                               IModel*                  model,
                               const DateTime&          eDate,
                               EventResults*            events) const = 0;
    };

    CClassConstSP getEventInterface() const;

    virtual void getEvents(IObjectConstSP       obj,
                           IModel*              model,
                           const DateTime&      eDate,
                           EventResults*        events) const;

    // if we get risk and legal events we need to know whether a risk
    // event will just be a duplicated legal event (and hence removable)
    virtual bool isValidRiskEvent(const DateTime& today) const;

private:
    KnownCashflows();
    KnownCashflows(const KnownCashflows &rhs);
    KnownCashflows& operator=(const KnownCashflows& rhs);
    static IObject* defaultKnownCashflows();
    static void load(CClassSP& clazz);

    CashFlowArraySP   cashFlows;
    string            ccy;
	CashflowInfoArraySP info; // optional
};
typedef smartPtr<KnownCashflows> KnownCashflowsSP;


/*************************/

class RISKMGR_DLL AllCashflows : public Event, public virtual IEvent{
public:
    static CClassConstSP const TYPE;

    AllCashflows(const DateTime& eDate,
        CashflowInfoArraySP cashflowInfos, string ccy);

    // interface instrument must satisfy to handle this event
    class RISKMGR_DLL IEventHandler {
    public:
        static CClassConstSP const TYPE;

        virtual void getEvents(const AllCashflows*    flows,
                               IModel*                  model,
                               const DateTime&          eDate,
                               EventResults*            events) const = 0;
    };

    CClassConstSP getEventInterface() const;

    virtual void getEvents(IObjectConstSP       obj,
                           IModel*              model,
                           const DateTime&      eDate,
                           EventResults*        events) const;

    // if we get risk and legal events we need to know whether a risk
    // event will just be a duplicated legal event (and hence removable)
    virtual bool isValidRiskEvent(const DateTime& today) const;

private:
    AllCashflows();
    AllCashflows(const AllCashflows &rhs);
    AllCashflows& operator=(const AllCashflows& rhs);
    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

	CashflowInfoArraySP cashflowInfos;
    string ccy;
};
typedef smartPtr<AllCashflows> AllCashflowsSP;

/******************************************************************************/
/** Event (implemented by intrument) to provide the QLib client with a list of 
 *  pas fixings needed (before-pricing pass).
 */
class RISKMGR_DLL FixingReqEvent : public Event, public virtual IEvent{
public:

    // interface instrument must satisfy to handle this event
    class RISKMGR_DLL IEventHandler {
    public:
        static CClassConstSP const TYPE;
        virtual void getEvents(const FixingReqEvent*,
                               IModel*             model,
                               const DateTime&     eDate,
                               EventResults*       events) const = 0;
    };

    virtual void getEvents(IObjectConstSP       obj,
                           IModel*              model,
                           const DateTime&      eDate,
                           EventResults*        events) const;


    // if we get risk and legal events we need to know whether a risk
    // event will just be a duplicated legal event (and hence removable)
    virtual bool isValidRiskEvent(const DateTime& today) const;

    CClassConstSP getEventInterface() const { return IEventHandler::TYPE; }
    FixingReqEvent(const CClassConstSP &type=TYPE) : Event(type) {}
    static CClassConstSP const TYPE;

    string indexSpecName;
    DateTimeArray dates;

private:
    FixingReqEvent(const FixingReqEvent &rhs);
    FixingReqEvent& operator=(const FixingReqEvent& rhs);
    static IObject* defaultConstructor(void) { return new FixingReqEvent(); }
    static void load(CClassSP& clazz);
};

typedef smartPtr<FixingReqEvent> FixingReqEventSP;
typedef array<FixingReqEventSP, FixingReqEvent> FixingReqEventArray;

DRLIB_END_NAMESPACE

#endif
