//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : CreditTree.cpp
//
//   Description : top level model class for rates trees
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CreditTree.hpp"

DRLIB_BEGIN_NAMESPACE 

CreditTree::CreditTree(CClassConstSP type): RateTree(type) {}

void CreditTree::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(CreditTree, clazz);
    SUPERCLASS(RateTree);

    // no fields yet
}

CreditTree::~CreditTree() {
}

FDProductSP CreditTree::makeProduct(const IProdCreatorSP & creator) {
    throw ModelException("Charles to implement ...", "CreditTree::makeProduct");
    return FDProductSP();
}

CClassConstSP const CreditTree::TYPE = CClass::registerClassLoadMethod(
    "CreditTree", typeid(CreditTree), CreditTree::load );

// for class loading
bool CreditTreeLoad(void) {
    return (CreditTree::TYPE != 0);
}

DRLIB_END_NAMESPACE
