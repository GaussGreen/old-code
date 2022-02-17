//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : NullDecretionCurve.cpp
//
//   Description : A null decretion curve which always returns principal balance as 1
//
//   Author      : Keith Jia
//
//   Date        : 03 January 2006
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/Class.hpp"
#include "edginc/NullDecretionCurve.hpp"
#include "edginc/RollingSettlement.hpp"

DRLIB_BEGIN_NAMESPACE

NullDecretionCurve::NullDecretionCurve(const string& inName)
    :DecretionCurve(TYPE)
{
    name = inName;
    prepaySettle = RollingSettlementSP(new RollingSettlement());
}

NullDecretionCurve::NullDecretionCurve(CClassConstSP clazz)
    :DecretionCurve(TYPE)
{}

NullDecretionCurve::~NullDecretionCurve()
{}

//------------------------------------------------
// Override IDecretionCurve pure vitual functions
//------------------------------------------------

/** pv, which actually is balance here, always 1 */
double NullDecretionCurve::pv(const DateTime& startDate,
                              const DateTime& endDate) const 
{
    return 1.0;
}
    
/** pv (balance) relative to initial balance */
double NullDecretionCurve::pv(const DateTime& endDate) const {
    return 1.0;
}

/** return decretion speed on a date */
double NullDecretionCurve::getDecretionSpeed(const DateTime& date) const {
    return 0.0;
}

/** return if balances are stepwise or continuous */ 
bool NullDecretionCurve::isStepBalances() const { 
    return true; 
}

/** Returns a reference to the value date */
double NullDecretionCurve::getFactor(const DateTime&) const {
    return 1.0;
}

DateTimeArraySP NullDecretionCurve::getStepDates() const {
    return DateTimeArraySP( new DateTimeArray(0)  );
}

DateTimeArraySP NullDecretionCurve::getAlternativeStepDates() const {
    return DateTimeArraySP( new DateTimeArray(0)  );
}

SettlementConstSP NullDecretionCurve::getSettlement() const
{
    return prepaySettle;
}

//--------------------------------------------
// load
//--------------------------------------------

void NullDecretionCurve::load(CClassSP& clazz){
    clazz->setPublic();
    REGISTER(NullDecretionCurve, clazz);
    SUPERCLASS(DecretionCurve);
    EMPTY_SHELL_METHOD(defaultNullDecretionCurve);
    FIELD(prepaySettle, "prepay settle rule, default to immeidiate settle");
    FIELD_MAKE_TRANSIENT(prepaySettle);
}

IObject* NullDecretionCurve::defaultNullDecretionCurve(){
    return new NullDecretionCurve();
}


//--------------------------------------------
// TYPE
//--------------------------------------------
CClassConstSP const NullDecretionCurve::TYPE = CClass::registerClassLoadMethod(
    "NullDecretionCurve", typeid(NullDecretionCurve), load);

DRLIB_END_NAMESPACE
