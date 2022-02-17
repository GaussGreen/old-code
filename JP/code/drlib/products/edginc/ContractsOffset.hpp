//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ContractsOffset.hpp
//
//   Description : Contract Offset used for fixing.
//
//   Author      : Simon A Creeger
//
//   Date        : 11 Dec 2006
//
//
//----------------------------------------------------------------------------

#ifndef EDG_CONTRACTSOFFSET_HPP
#define EDG_CONTRACTSOFFSET_HPP

#include "edginc/Object.hpp"
#include "edginc/Array.hpp"
#include "edginc/AtomicArray.hpp"

DRLIB_BEGIN_NAMESPACE

class ContractsOffset;
typedef smartConstPtr<ContractsOffset> ContractsOffsetConstSP;
typedef smartPtr<ContractsOffset> ContractsOffsetSP;

typedef array<ContractsOffsetSP,ContractsOffset> ContractsOffsetArray;
typedef smartPtr<ContractsOffsetArray> ContractsOffsetArraySP;
typedef smartConstPtr<ContractsOffsetArray> ContractsOffsetArrayConstSP;

typedef array<ContractsOffsetArraySP,ContractsOffsetArray> ContractsOffsetArrayArray;
typedef smartPtr<ContractsOffsetArrayArray> ContractsOffsetArrayArraySP;
typedef smartConstPtr<ContractsOffsetArrayArray> ContractsOffsetArrayArrayConstSP;

class PRODUCTS_DLL ContractsOffset : public CObject {
public:
    static CClassConstSP const TYPE;

    struct offsetAssetType {
        enum Enum {ENERGY, METAL, OTHER};
    };
    typedef BoxedEnum<offsetAssetType::Enum> offsetAssetTypeBoxedEnum;

    offsetAssetType::Enum getAsset() const;

    int getOffset() const;

    static void ContractsOffset::load(CClassSP& clazz) {
        REGISTER(ContractsOffset, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultContractsOffset);
        FIELD(offsetAsset, "Enum Asset [ENERGY(1),METAL(2),OTHER(3)]");
        FIELD(offset, "Contract offset index");  
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

private:
    ContractsOffset();

    offsetAssetType::Enum      offsetAsset;
    int                        offset;

    static IObject* defaultContractsOffset();
};

DRLIB_END_NAMESPACE
#endif
