/**
 * @file CreditIndexBasisParallel.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/RiskProperty_TYPES.hpp"
#include "edginc/CreditIndexBasisParallel.hpp"

DRLIB_BEGIN_NAMESPACE

CreditIndexBasisParallel::CreditIndexBasisParallel(): CObject(TYPE) {}
CreditIndexBasisParallel::~CreditIndexBasisParallel() {}

static void CreditIndexBasisParallel_load(CClassSP& clazz) {
    REGISTER(CreditIndexBasisParallel, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(DefaultConstructor<CreditIndexBasisParallel>::iObject);
}

CClassConstSP const CreditIndexBasisParallel::TYPE = CClass::registerClassLoadMethod("CreditIndexBasisParallel", typeid(CreditIndexBasisParallel), CreditIndexBasisParallel_load);

RiskProperty_TYPES(CreditIndexBasisParallel)

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
