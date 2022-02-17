//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : PrepayParallelTP.cpp
//
//   Description : ABCDS prepay curve horizontal parallel shift Tweak Property
//
//   Date        : May 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/RiskProperty_TYPES.hpp"
#include "edginc/PrepayParallelTP.hpp"

DRLIB_BEGIN_NAMESPACE

PrepayParallelTP::PrepayParallelTP(): CObject(TYPE) {}
PrepayParallelTP::~PrepayParallelTP() {}

static void PrepayParallelTP_load(CClassSP& clazz) {
    REGISTER(PrepayParallelTP, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(DefaultConstructor<PrepayParallelTP>::iObject);
}

CClassConstSP const PrepayParallelTP::TYPE = CClass::registerClassLoadMethod("PrepayParallelTP", typeid(PrepayParallelTP), PrepayParallelTP_load);

RiskProperty_TYPES(PrepayParallelTP)

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
