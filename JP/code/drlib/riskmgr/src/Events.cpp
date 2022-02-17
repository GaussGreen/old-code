#include "edginc/config.hpp"
#include "edginc/Events.hpp"
#include "edginc/Format.hpp"

DRLIB_BEGIN_NAMESPACE

// for reflection 
Event::Event(CClassConstSP clazz):  CObject(clazz) {}

// for inheritance 
Event::Event(CClassConstSP clazz, const DateTime& eDate): 
            CObject(clazz), eventDate(eDate) {}

// Invoked when Class is 'loaded'
void Event::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(Event, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultEvent);
    FIELD(eventDate, "The date the event occurred");
//    FIELD(eventAction, "Any actions required by this event");
}

IObject* Event::defaultEvent(){
    return new Event(TYPE);
}

CClassConstSP const Event::TYPE = CClass::registerClassLoadMethod(
    "Event", typeid(Event), load);

DEFINE_TEMPLATE_TYPE(EventArray);

CClassConstSP const IEvent::TYPE = CClass::registerInterfaceLoadMethod(
    "IEvent", typeid(IEvent), 0);

/******************************************************************************/

const string BarrierBreach::KNOCK_IN = "Knock In";
const string BarrierBreach::KNOCK_OUT = "Knock Out";
const string BarrierBreach::LOCK_IN = "Lock In";
const string BarrierBreach::NOT_APPLICABLE = "Not Applicable";
const string BarrierBreach::DAILY = "Daily";
const string BarrierBreach::EUROPEAN = "European";
const string BarrierBreach::CONTINUOUS = "Continuous";

/* for reflection */
BarrierBreach::BarrierBreach(): Event(TYPE) {}

// Invoked when Class is 'loaded'
void BarrierBreach::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(BarrierBreach, clazz);
    SUPERCLASS(Event);
    IMPLEMENTS(IEvent);
    EMPTY_SHELL_METHOD(defaultBarrierBreach);
    FIELD(barrDesc, "Description of barrier breached");
    FIELD(barrierType, "The type of barrier breached");
    FIELD(monType, "The monitoring type of the barrier - european/daily/cts");
    FIELD(isUp, "Is the barrier up?");
    FIELD(hitsRemaining, "How many hits still to come?");
    FIELD(names, "The assets/quantities that breached");
    FIELD(levels, "The levels of the assets/quantities that breached");
    FIELD(barrierLevels, "The barrier levels for the assets/quantities that breached");
}

BarrierBreach::BarrierBreach(const DateTime& eDate, string barrDesc, 
                             string monType, string barrType, bool isItUp, 
                             int hitsLeft, StringArraySP names, 
                             DoubleArraySP levels, DoubleArraySP bLevels) :
        Event(TYPE, eDate), barrDesc(barrDesc), monType(monType), 
        barrierType(barrType), isUp(isItUp),
        hitsRemaining(hitsLeft), names(names), 
        levels(levels), barrierLevels(bLevels) {};

IObject* BarrierBreach::defaultBarrierBreach(){
    return new BarrierBreach();
}

CClassConstSP BarrierBreach::getEventInterface() const {
    return IEventHandler::TYPE;
}

void BarrierBreach::getEvents(IObjectConstSP       obj,
                              IModel*              model,
                              const DateTime&      eDate,
                              EventResults*        events) const {
    const IEventHandler& eventObj = dynamic_cast<const IEventHandler&>(*obj);
    eventObj.getEvents(this, model, eDate, events);
}

// if we get risk and legal events we need to know whether a risk
// event will just be a duplicated legal event (and hence removable)
bool BarrierBreach::isValidRiskEvent(const DateTime& today) const {
    // risk barrier breaches strictly before today would be duplicated as legal
    return (today <= eventDate);
}

CClassConstSP const BarrierBreach::TYPE = CClass::registerClassLoadMethod(
    "BarrierBreach", typeid(BarrierBreach), load);

CClassConstSP const BarrierBreach::IEventHandler::TYPE = CClass::registerInterfaceLoadMethod(
    "BarrierBreach::IEventHandler", typeid(BarrierBreach::IEventHandler), 0);

/******************************************************************************/

/* for reflection */
FlexBarrierBreach::FlexBarrierBreach(): Event(TYPE) {}

// Invoked when Class is 'loaded'
void FlexBarrierBreach::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(FlexBarrierBreach, clazz);
    SUPERCLASS(Event);
    IMPLEMENTS(IEvent);
    EMPTY_SHELL_METHOD(defaultFlexBarrierBreach);
    FIELD(variableName, "The name of the variable being monitored");
    FIELD(isUp, "Is the barrier up?");
    FIELD(variableLevel, "The value of the variable at breach");
    FIELD(barrierLevel, "The value of the barrier at breach");
}

FlexBarrierBreach::FlexBarrierBreach(const DateTime& eDate,
                                     string varName, double barrLevel,
                                     double varLevel, bool isItUp) :
     Event(TYPE, eDate), variableName(varName), isUp(isItUp), variableLevel(varLevel),
     barrierLevel(barrLevel){};

IObject* FlexBarrierBreach::defaultFlexBarrierBreach(){
    return new FlexBarrierBreach();
}

CClassConstSP FlexBarrierBreach::getEventInterface() const {
    return IEventHandler::TYPE;
}

void FlexBarrierBreach::getEvents(IObjectConstSP       obj,
                                  IModel*              model,
                                  const DateTime&      eDate,
                                  EventResults*        events) const {
    const IEventHandler& eventObj = dynamic_cast<const IEventHandler&>(*obj);
    eventObj.getEvents(this, model, eDate, events);
}

// if we get risk and legal events we need to know whether a risk
// event will just be a duplicated legal event (and hence removable)
bool FlexBarrierBreach::isValidRiskEvent(const DateTime& today) const {
    // risk barrier breaches strictly before today would be duplicated as legal
    return (today <= eventDate);
}

CClassConstSP const FlexBarrierBreach::TYPE = CClass::registerClassLoadMethod(
    "FlexBarrierBreach", typeid(FlexBarrierBreach), load);

CClassConstSP const FlexBarrierBreach::IEventHandler::TYPE = CClass::registerInterfaceLoadMethod(
    "FlexBarrierBreach::IEventHandler", typeid(FlexBarrierBreach::IEventHandler), 0);

/******************************************************************************/

const string TargetRedemption::KNOCK_IN = "Knock In";
const string TargetRedemption::KNOCK_OUT = "Knock Out";
const string TargetRedemption::NOT_APPLICABLE = "Not Applicable";

/* for reflection */
TargetRedemption::TargetRedemption(): Event(TYPE) {}

// Invoked when Class is 'loaded'
void TargetRedemption::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(TargetRedemption, clazz);
    SUPERCLASS(Event);
    IMPLEMENTS(IEvent);
    EMPTY_SHELL_METHOD(defaultTargetRedemption);
    FIELD(coupon, "The breaching coupon");
    FIELD(sumCoupons, "The sum of the coupons so far");
    FIELD(targetLevel, "The target level");
    FIELD(bonus, "The bonus coupon");
    FIELD(redemption, "The total redemption value");
    FIELD(floatingLeg, "The floating leg behaviour");
}

TargetRedemption::TargetRedemption(const DateTime& eDate,
                                   double cpn, double sumCpn,
                                   double target, double bonusCpn,
                                   double redemptionVal, string floatingLegType) :
     Event(TYPE, eDate), coupon(cpn), sumCoupons(sumCpn), 
     targetLevel(target), bonus(bonusCpn), redemption(redemptionVal),
     floatingLeg(floatingLegType){};

IObject* TargetRedemption::defaultTargetRedemption(){
    return new TargetRedemption();
}

CClassConstSP TargetRedemption::getEventInterface() const {
    return IEventHandler::TYPE;
}

void TargetRedemption::getEvents(IObjectConstSP       obj,
                                 IModel*              model,
                                 const DateTime&      eDate,
                                 EventResults*        events) const {
    const IEventHandler& eventObj = dynamic_cast<const IEventHandler&>(*obj);
    eventObj.getEvents(this, model, eDate, events);
}

// if we get risk and legal events we need to know whether a risk
// event will just be a duplicated legal event (and hence removable)
bool TargetRedemption::isValidRiskEvent(const DateTime& today) const {
    // risk early redemptions strictly before today would be duplicated as legal
    return (today <= eventDate);
}

CClassConstSP const TargetRedemption::TYPE = CClass::registerClassLoadMethod(
    "TargetRedemption", typeid(TargetRedemption), load);

CClassConstSP const TargetRedemption::IEventHandler::TYPE = CClass::registerInterfaceLoadMethod(
    "TargetRedemption::IEventHandler", typeid(TargetRedemption::IEventHandler), 0);

/******************************************************************************/

const string Callability::CALLABLE = "Callable";
const string Callability::PUTTABLE = "Puttable";
const string Callability::CALL     = "Call";
const string Callability::PUT      = "Put";

/* for reflection */
Callability::Callability(): Event(TYPE) {}

// Invoked when Class is 'loaded'
void Callability::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(Callability, clazz);
    SUPERCLASS(Event);
    IMPLEMENTS(IEvent);
    EMPTY_SHELL_METHOD(defaultCallability);
    FIELD(notifDate, "The date on which the decision must be made");
    FIELD(callDate, "The date on which the call/put is exercised");
    FIELD(type, "Callable or puttable");
    FIELD(amount, "The amount due in order to exercise");
}

Callability::Callability(const DateTime& notifDate, const DateTime& callDate, 
                         string type, double amount) :
     Event(TYPE, notifDate), notifDate(notifDate), callDate(callDate), 
     type(type), amount(amount) {};

IObject* Callability::defaultCallability(){
    return new Callability();
}

CClassConstSP Callability::getEventInterface() const {
    return IEventHandler::TYPE;
}

void Callability::getEvents(IObjectConstSP       obj,
                            IModel*              model,
                            const DateTime&      eDate,
                            EventResults*        events) const {
    const IEventHandler& eventObj = dynamic_cast<const IEventHandler&>(*obj);
    eventObj.getEvents(this, model, eDate, events);
}

// if we get risk and legal events we need to know whether a risk
// event will just be a duplicated legal event (and hence removable)
bool Callability::isValidRiskEvent(const DateTime& today) const {
    // risk callability is always a duplicate of legal
    return false;
}

CClassConstSP const Callability::TYPE = CClass::registerClassLoadMethod(
    "Callability", typeid(Callability), load);

CClassConstSP const Callability::IEventHandler::TYPE = CClass::registerInterfaceLoadMethod(
    "Callability::IEventHandler", typeid(Callability::IEventHandler), 0);

/******************************************************************************/

const string SPIFixing::ZC_BOND_FIXING = "ZC Bond Fixing";
const string SPIFixing::MM_RATE_FIXING = "Loan Cost Fixing";
const string SPIFixing::BOND_FLOOR_FIXING = "Bond Floor Fixing";

/* for reflection */
SPIFixing::SPIFixing(): Event(TYPE) {}

// Invoked when Class is 'loaded'
void SPIFixing::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(SPIFixing, clazz);
    SUPERCLASS(Event);
    IMPLEMENTS(IEvent);
    EMPTY_SHELL_METHOD(defaultSPIFixing);
    FIELD(type, "ZC Bond Fixing, Loan Cost Fixing or Bond Floor Fixing");
    FIELD(amount, "The fixing");
}

SPIFixing::SPIFixing(const DateTime& eDate, string type, double amount) :
     Event(TYPE, eDate), type(type), amount(amount) {};

IObject* SPIFixing::defaultSPIFixing(){
    return new SPIFixing();
}

CClassConstSP SPIFixing::getEventInterface() const {
    return IEventHandler::TYPE;
}

void SPIFixing::getEvents(IObjectConstSP       obj,
                          IModel*              model,
                          const DateTime&      eDate,
                          EventResults*        events) const {
    const IEventHandler& eventObj = dynamic_cast<const IEventHandler&>(*obj);
    eventObj.getEvents(this, model, eDate, events);
}

// if we get risk and legal events we need to know whether a risk
// event will just be a duplicated legal event (and hence removable)
bool SPIFixing::isValidRiskEvent(const DateTime& today) const {
    // no concept of risk/legal for SPIFixings
    return false;
}

CClassConstSP const SPIFixing::TYPE = CClass::registerClassLoadMethod(
    "SPIFixing", typeid(SPIFixing), load);

CClassConstSP const SPIFixing::IEventHandler::TYPE = CClass::registerInterfaceLoadMethod(
    "SPIFixing::IEventHandler", typeid(SPIFixing::IEventHandler), 0);


/******************************************************************************/

void CashflowInfo::merge(const CashflowInfo &cfi, string const &prefix) {
    if (cfi.amountType == CfType::UNSET) {
        // nothing new to add
        return;
    }
    updateAmountType(cfi.amountType);
    if (amountType == AmountType::UNKNOWN) {
        // do not add key/values, it is irrelevant anyway
        return;
    }
    if (prefix.size()) {
        for (int i=0; i < cfi.keys.size(); ++i) {
            keys.push_back(prefix + "." + cfi.keys[i]);
        }
    }
    else keys.insert(keys.end(), cfi.keys.begin(), cfi.keys.end());
    values.insert(values.end(), cfi.values.begin(), cfi.values.end());
}

void CashflowInfo::push(const string &key, const string &value) { 
	keys.push_back(key); 
	values.push_back(value); 
}

void CashflowInfo::push(const string &key, double value) {
	push(key, Format::toString(value)); 
}

void CashflowInfo::push(const string &key, bool value) {
	push(key, Format::toString(value)); 
}

void CashflowInfo::push(const string &key, DateTime value) {
	push(key, value.toString());
}

void CashflowInfo::setMultiplier(double multiplier) {
    amount *= multiplier;
    multiplierUsed *= multiplier;
}

CashflowInfo::CashflowInfo(const CClassConstSP &type) : CObject(type), 
cfType(CfType::UNSET), amount(0.), amountType(AmountType::UNSET), multiplierUsed(1.) {}

IObject* CashflowInfo::defaultCashflowInfo(void) {
	return new CashflowInfo();
}


void CashflowInfo::updateAmountType(AmountType::Enum newAmountType) {
    if (amountType == AmountType::UNSET) {
        amountType = newAmountType;
        return;
    }
    if (amountType == AmountType::UNKNOWN) {
        // UNKNOWN cannot be updgraded to anything else
        return;
    }
    switch (newAmountType) {
        case AmountType::UNSET:
            throw ModelException(__FUNCTION__, "UNSET: unauthorized value");
        case AmountType::KNOWN:
            break;
        case AmountType::ESTIMATED:
            if (amountType == AmountType::KNOWN)
                amountType = AmountType::ESTIMATED;
            break;
        case AmountType::UNKNOWN:
            amountType = AmountType::UNKNOWN;
            keys.clear();
            values.clear();
            amount = 0.;
            break;
        default: ASSERT(0);
    }
}

// Invoked when Class is 'loaded'
void CashflowInfo::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(CashflowInfo, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultCashflowInfo);
    FIELD(componentName, "");
    FIELD(cfType, "");
    FIELD(amount,"");
    FIELD(date,"");
    FIELD(amountType, "");
    FIELD(multiplierUsed, "");
    FIELD(keys, "keys");
    FIELD(values, "values");
}

CClassConstSP const CashflowInfo::TYPE = CClass::registerClassLoadMethod(
    "CashflowInfo", typeid(CashflowInfo), CashflowInfo::load);

DEFINE_TEMPLATE_TYPE(CashflowInfoArray);


START_PUBLIC_ENUM_DEFINITION(CashflowInfo::CfType::Enum, "Cashflow type");
ENUM_VALUE_AND_NAME(CashflowInfo::CfType::UNSET, "UNSET", "Unknown type");
ENUM_VALUE_AND_NAME(CashflowInfo::CfType::PRINCIPAL, "PRINCIPAL", "Principal payment");
ENUM_VALUE_AND_NAME(CashflowInfo::CfType::COUPON, "COUPON", "Coupon payment");
ENUM_VALUE_AND_NAME(CashflowInfo::CfType::REDEEMER, "REDEEMER", "Redeemer payment");
END_ENUM_DEFINITION(CashflowInfo::CfType::Enum);

START_PUBLIC_ENUM_DEFINITION(CashflowInfo::AmountType::Enum, "Cashflow type");
ENUM_VALUE_AND_NAME(CashflowInfo::AmountType::UNSET, "UNSET", "Unknown type");
ENUM_VALUE_AND_NAME(CashflowInfo::AmountType::KNOWN, "KNOWN", "Amount known for sure");
ENUM_VALUE_AND_NAME(CashflowInfo::AmountType::ESTIMATED, "ESTIMATED", "Amount estimated (rough approximation)");
ENUM_VALUE_AND_NAME(CashflowInfo::AmountType::UNKNOWN, "UNKNOWN", "A payment will occur at this date but the amount is not known yet");
END_ENUM_DEFINITION(CashflowInfo::AmountType::Enum);

/******************************************************************************/
/* for reflection */
KnownCashflows::KnownCashflows(): Event(TYPE) {}

// Invoked when Class is 'loaded'
void KnownCashflows::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(KnownCashflows, clazz);
    SUPERCLASS(Event);
    IMPLEMENTS(IEvent);
    EMPTY_SHELL_METHOD(defaultKnownCashflows);
    FIELD(cashFlows, "cashflows");
    FIELD(info, "info"); FIELD_MAKE_OPTIONAL(info);
    FIELD(ccy, "currency");
}

KnownCashflows::KnownCashflows(const DateTime& eDate, CashFlowArraySP cashFlows, 
                               string ccy, CashflowInfoArraySP info) :
     Event(TYPE, eDate), cashFlows(cashFlows), ccy(ccy), info(info) {};

IObject* KnownCashflows::defaultKnownCashflows(){
    return new KnownCashflows();
}

CClassConstSP KnownCashflows::getEventInterface() const {
    return IEventHandler::TYPE;
}

void KnownCashflows::getEvents(IObjectConstSP       obj,
                               IModel*              model,
                               const DateTime&      eDate,
                               EventResults*        events) const {
    const IEventHandler& eventObj = dynamic_cast<const IEventHandler&>(*obj);
    eventObj.getEvents(this, model, eDate, events);
}

// if we get risk and legal events we need to know whether a risk
// event will just be a duplicated legal event (and hence removable)
bool KnownCashflows::isValidRiskEvent(const DateTime& today) const {
    // no concept of risk/legal for known cashflows at the moment
    // long term maybe delegate to the context
    // for example a rebate triggered by a barrier breach might be a risk event
    // but a fixed flow from a Libor leg wouldn't be
    return false;
}

CClassConstSP const KnownCashflows::TYPE = CClass::registerClassLoadMethod(
    "KnownCashflows", typeid(KnownCashflows), load);

CClassConstSP const KnownCashflows::IEventHandler::TYPE = CClass::registerInterfaceLoadMethod(
    "KnownCashflows::IEventHandler", typeid(KnownCashflows::IEventHandler), 0);

/******************************************************************************/

AllCashflows::AllCashflows() : Event(TYPE) {}

AllCashflows::AllCashflows(
    const DateTime& eDate, CashflowInfoArraySP cashflowInfos, string ccy) 
: Event(TYPE, eDate), cashflowInfos(cashflowInfos), ccy(ccy) 
{}

// Invoked when Class is 'loaded'
void AllCashflows::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(AllCashflows, clazz);
    SUPERCLASS(Event);
    IMPLEMENTS(IEvent);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(cashflowInfos, "cashflowInfos");
    FIELD(ccy, "currency");
}


IObject* AllCashflows::defaultConstructor(){
    return new AllCashflows();
}

CClassConstSP AllCashflows::getEventInterface() const {
    return IEventHandler::TYPE;
}

void AllCashflows::getEvents(IObjectConstSP       obj,
                               IModel*              model,
                               const DateTime&      eDate,
                               EventResults*        events) const {
    const IEventHandler& eventObj = dynamic_cast<const IEventHandler&>(*obj);
    eventObj.getEvents(this, model, eDate, events);
}

// if we get risk and legal events we need to know whether a risk
// event will just be a duplicated legal event (and hence removable)
bool AllCashflows::isValidRiskEvent(const DateTime& today) const {
    // no concept of risk/legal for known cashflows at the moment
    // long term maybe delegate to the context
    // for example a rebate triggered by a barrier breach might be a risk event
    // but a fixed flow from a Libor leg wouldn't be
    return false;
}

CClassConstSP const AllCashflows::TYPE = CClass::registerClassLoadMethod(
    "AllCashflows", typeid(AllCashflows), load);

CClassConstSP const AllCashflows::IEventHandler::TYPE = CClass::registerInterfaceLoadMethod(
    "AllCashflows::IEventHandler", typeid(AllCashflows::IEventHandler), 0);


/******************************** FixingReqEvent *****************************/
// Invoked when Class is 'loaded'
void FixingReqEvent::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(FixingReqEvent, clazz);
    SUPERCLASS(Event);
    IMPLEMENTS(IEvent);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(indexSpecName, "");
    FIELD(dates, "");
}

void FixingReqEvent::getEvents(IObjectConstSP      obj,
                              IModel*              model,
                              const DateTime&      eDate,
                              EventResults*        events
) const {
    const IEventHandler& eventObj = dynamic_cast<const IEventHandler&>(*obj);
    eventObj.getEvents(this, model, eDate, events);
}

// if we get risk and legal events we need to know whether a risk
// event will just be a duplicated legal event (and hence removable)
bool FixingReqEvent::isValidRiskEvent(const DateTime& today) const {
    // no concept of risk/legal for Fixing events
    return false;
}

CClassConstSP const FixingReqEvent::TYPE = CClass::registerClassLoadMethod(
    "FixingReqEvent", typeid(FixingReqEvent), FixingReqEvent::load);

CClassConstSP const FixingReqEvent::IEventHandler::TYPE = CClass::registerInterfaceLoadMethod(
    "FixingReqEvent::IEventHandler", typeid(FixingReqEvent::IEventHandler), 0);

DRLIB_END_NAMESPACE
