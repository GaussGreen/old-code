//----------------------------------------------------------------------------
//
//   Group       : New York Credit QR
//
//   Filename    : KRiskyComponent.cpp
//
//   Description : generic instrument for components with credit risk
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/ICDSParSpreads.hpp"
#include "edginc/KRiskyComponent.hpp"

DRLIB_BEGIN_NAMESPACE

void KRiskyComponent::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    clazz->setDescription("QLIB product representation of risky instrument");
    REGISTER(KRiskyComponent, clazz);
    SUPERCLASS(KComponent);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(cdsParSpreads,"Credit curve");
    FIELD_MAKE_OPTIONAL(cdsParSpreads);
}

CClassConstSP const KRiskyComponent::TYPE = CClass::registerClassLoadMethod(
    "KRiskyComponent", typeid(KRiskyComponent), KRiskyComponent::load);

bool KRiskyComponentLoad(){
    return (KRiskyComponent::TYPE != 0);
}

DRLIB_END_NAMESPACE

