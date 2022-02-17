//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PhysicalAndCashSettlement.cpp
//
//   Description : special settlement designed for DailyAccumulator.cpp
//
//   Author      : Keith Law
//
//   Date        : 12 Dec 2006
//
//
//   $Log:
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/PhysicalAndCashSettlement.hpp"

DRLIB_BEGIN_NAMESPACE

PhysicalAndCashSettlement::PhysicalAndCashSettlement(): CObject(TYPE),
 isPhysicalAndCash(false),
 phyDeliveredAtSpot(false),
 physicalWeight(1.0),
 cashWeight(0.0){};


class PhysicalAndCashSettlementHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(PhysicalAndCashSettlement, clazz);
        //SUPERCLASS(Settlement);
        EMPTY_SHELL_METHOD(defaultPhysicalAndCashSettlement);
        FIELD(isPhysicalAndCash, "For physical and cash settlement purpose only");  
        FIELD(physicalWeight, "Weight for physical settlement");  
        FIELD_MAKE_OPTIONAL(physicalWeight);
        FIELD(cashWeight, "Weight for cash settlement");  
        FIELD_MAKE_OPTIONAL(cashWeight);
        FIELD(phyDeliveredAtSpot, "if Physical settle, delivered at spot price?");  
    }

    static IObject* defaultPhysicalAndCashSettlement(){
        return new PhysicalAndCashSettlement();
    }
};



CClassConstSP const PhysicalAndCashSettlement::TYPE = CClass::registerClassLoadMethod(
    "PhysicalAndCashSettlement", typeid(PhysicalAndCashSettlement), PhysicalAndCashSettlementHelper::load);

DRLIB_END_NAMESPACE
