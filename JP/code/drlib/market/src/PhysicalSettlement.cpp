//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PhysicalSettlement.cpp
//
//   Description : Physical style instrument settlement
//
//   Author      : Andrew J Swain
//
//   Date        : 27 February 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/PhysicalSettlement.hpp"
#include "edginc/Addin.hpp"


DRLIB_BEGIN_NAMESPACE

PhysicalSettlement::PhysicalSettlement(): InstrumentSettlement(TYPE) {}

PhysicalSettlement::~PhysicalSettlement() {}

/** given a trade date, when does this settle ? */
DateTime PhysicalSettlement::settles(const DateTime& date,
                                     const CAsset*   asset) const {
static const string method = "PhysicalSettlement::settles";
try {
    if (!asset) {
        throw ModelException(method, "Unable to compute a physical settlement date"
                                     " when no asset has been passed");
    }

    if (!asset->canPhysicallySettle()) {
        throw ModelException("PhysicalSettlement::settles",
                             "Physical settlement is disallowed for assets "
                             "of type " + asset->getClass()->getName() + 
                             " (asset is " + asset->getName() + ")");
    }

    return asset->settleDate(date);
}
catch (exception& e) {
    throw ModelException(e, method);
}
}

/** is this a physical settlement ? */
bool PhysicalSettlement::isPhysical() const {
    return true;
}

/** is this a margin option ? */
bool PhysicalSettlement::isMargin() const {
    return false;
}

class PhysicalSettlementHelper{
public:
    /** Invoked when PhysicalSettlement is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(PhysicalSettlement, clazz);
        SUPERCLASS(InstrumentSettlement);
        EMPTY_SHELL_METHOD(defaultPhysicalSettlement);
        Addin::registerConstructor(Addin::UTILITIES, PhysicalSettlement::TYPE);
    }

    static IObject* defaultPhysicalSettlement(){
        return new PhysicalSettlement();
    }
};

CClassConstSP const PhysicalSettlement::TYPE = CClass::registerClassLoadMethod(
    "PhysicalSettlement", typeid(PhysicalSettlement), PhysicalSettlementHelper::load);

DRLIB_END_NAMESPACE
