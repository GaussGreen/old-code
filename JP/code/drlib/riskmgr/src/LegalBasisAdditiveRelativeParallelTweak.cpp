/**
 * @file LegalBasisAdditiveRelativeParallelTweak.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#define QLIB_LegalBasisAdditiveRelativeParallelTweak_CPP
#include "edginc/Expiry.hpp"
#include "edginc/RiskProperty_TYPES.hpp"
#include "edginc/LegalBasisAdditiveRelativeParallelTweak.hpp"

DRLIB_BEGIN_NAMESPACE

LegalBasisAdditiveRelativeParallelTweak::LegalBasisAdditiveRelativeParallelTweak(): CObject(TYPE) {}
LegalBasisAdditiveRelativeParallelTweak::~LegalBasisAdditiveRelativeParallelTweak() {}

static void LegalBasisAdditiveRelativeParallelTweak_load(CClassSP& clazz) {
    REGISTER(LegalBasisAdditiveRelativeParallelTweak, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(DefaultConstructor<LegalBasisAdditiveRelativeParallelTweak>::iObject);
}

CClassConstSP const LegalBasisAdditiveRelativeParallelTweak::TYPE = CClass::registerClassLoadMethod("LegalBasisAdditiveRelativeParallelTweak", typeid(LegalBasisAdditiveRelativeParallelTweak), LegalBasisAdditiveRelativeParallelTweak_load);

RiskProperty_TYPES(LegalBasisAdditiveRelativeParallelTweak)

INSTANTIATE_TEMPLATE(class RISKMGR_DLL RiskProperty<LegalBasisAdditiveRelativeParallelTweak>);

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
