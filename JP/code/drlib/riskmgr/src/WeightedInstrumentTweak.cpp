//----------------------------------------------------------------------------
//
//   Group       : QR - Credit Hybrids
//
//   Filename    : WeightedInstrumentTweak.cpp
//
//   Description : Tag class for shift in instrument weight 
//                 (WeightedInstrumentCollection)
//
//   Author      : Linus Thand
//
//   Date        : 13 July 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/RiskProperty_TYPES.hpp"
#include "edginc/WeightedInstrumentTweak.hpp"

DRLIB_BEGIN_NAMESPACE

WeightedInstrumentTweak::WeightedInstrumentTweak(IntArrayConstSP instruments):
    CObject(TYPE),
    instruments(instruments)
{}

WeightedInstrumentTweak::~WeightedInstrumentTweak() {}

static void WeightedInstrumentTweak_load(CClassSP& clazz) {
    REGISTER(WeightedInstrumentTweak, clazz);
    SUPERCLASS(CObject);
    FIELD(instruments, "instruments");
    EMPTY_SHELL_METHOD(DefaultConstructor<WeightedInstrumentTweak>::iObject);
}

CClassConstSP const WeightedInstrumentTweak::TYPE = CClass::registerClassLoadMethod("WeightedInstrumentTweak", typeid(WeightedInstrumentTweak), WeightedInstrumentTweak_load);

RiskProperty_TYPES(WeightedInstrumentTweak)

DRLIB_END_NAMESPACE
