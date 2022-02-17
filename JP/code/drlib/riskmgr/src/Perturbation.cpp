//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Perturbation.cpp
//
//   Description : Defines the ability to shift a named piece of data
//
//   Author      : Mark A Robson/Andrew J Swain
//
//   Date        : 19 Jan 2005
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/Perturbation.hpp"
#include "edginc/SensMgr.hpp"

DRLIB_BEGIN_NAMESPACE

/* Define some strings that can be used by multiple scenarios to configure
 * the type of tweak to be performed.
 * Added here to enforce a consistent use of theses strings, and to avoid
 * creating them on the fly wherever needed. */
const string IPerturbation::ABSOLUTE = "A"; // ie, value += shift
const string IPerturbation::RELATIVE = "R"; // ie, value *= (1 + shift)
const string IPerturbation::SET      = "S"; // ie, value  = shift

IPerturbation::IPerturbation(){}
IPerturbation::~IPerturbation(){}

void IPerturbation::load(CClassSP& clazz){
    REGISTER_INTERFACE(IPerturbation, clazz);
    clazz->setPublic();
    EXTENDS(IObject);
}

CClassConstSP const IPerturbation::TYPE = CClass::registerInterfaceLoadMethod(
    "IPerturbation", typeid(IPerturbation), load);

Perturbation::~Perturbation(){}

/** Returns this */
ITweakNameResolver* Perturbation::nameResolver(){
    return this;
}

/** returns the name identifying the market data to be shifted. Returns
    null if not set */
OutputNameConstSP Perturbation::getMarketDataName() const{
    return marketDataName;
}

/** Does nothing */
void Perturbation::reset(){}

/** apply this Perturbation to the object with the specified name contained
    within the supplied object */
bool Perturbation::findAndShift(IObjectSP         objectToShift, 
                                OutputNameConstSP name){
    marketDataName = name;
    SensMgr sensMgr(objectToShift);
    sensMgr.shift(this);
    marketDataName.reset();
    return sensMgr.getShiftStatus();
}

bool IPerturbation::applyBeforeGetMarket() const {
    return false;
}

Perturbation::Perturbation(CClassConstSP clazz): CObject(clazz){}

void Perturbation::load(CClassSP& clazz){
    clazz->setPublic();
    REGISTER(Perturbation, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IPerturbation);
    FIELD_NO_DESC(marketDataName);
    FIELD_MAKE_TRANSIENT(marketDataName);
}

CClassConstSP const Perturbation::TYPE = CClass::registerClassLoadMethod(
    "Perturbation", typeid(Perturbation), load);

DRLIB_END_NAMESPACE
