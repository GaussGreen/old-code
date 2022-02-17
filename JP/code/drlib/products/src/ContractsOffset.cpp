//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ContractsOffset.cpp
//
//   Description : Contract Offset used for fixing.
//
//   Author      : Simon A Creeger
//
//   Date        : 11 Dec 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_CONTRACTSOFFSET_CPP
#include "edginc/ContractsOffset.hpp"

DRLIB_BEGIN_NAMESPACE

ContractsOffset::ContractsOffset() : CObject(TYPE), offsetAsset(ContractsOffset::offsetAssetType::OTHER), offset(0) {}

ContractsOffset::offsetAssetType::Enum ContractsOffset::getAsset() const {
    return offsetAsset;
}

int ContractsOffset::getOffset() const {
    return offset;
}

IObject* ContractsOffset::defaultContractsOffset() {
    return new ContractsOffset();
}

CClassConstSP const ContractsOffset::TYPE = CClass::registerClassLoadMethod(
    "ContractsOffset", typeid(ContractsOffset), ContractsOffset::load);

START_PUBLIC_ENUM_DEFINITION(ContractsOffset::offsetAssetType::Enum,"");
ENUM_VALUE_AND_NAME(ContractsOffset::offsetAssetType::ENERGY, "ENERGY","");
ENUM_VALUE_AND_NAME(ContractsOffset::offsetAssetType::METAL, "METAL","");
ENUM_VALUE_AND_NAME(ContractsOffset::offsetAssetType::OTHER, "OTHER","");
END_ENUM_DEFINITION(ContractsOffset::offsetAssetType::Enum);

DEFINE_TEMPLATE_TYPE(ContractsOffsetArray);
DEFINE_TEMPLATE_TYPE(ContractsOffsetArrayArray);

DRLIB_END_NAMESPACE