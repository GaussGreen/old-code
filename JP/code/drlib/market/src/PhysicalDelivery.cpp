//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PhysicalDelivery.cpp
//
//   Description : Results data holder for what needs to be physically settled
//
//   Author      : Andrew J Swain
//
//   Date        : 26 February 2004
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/PhysicalDelivery.hpp"
#include "edginc/AssetUtil.hpp"

DRLIB_BEGIN_NAMESPACE

// when physical should drop out of instrument price
DateTime PhysicalDelivery::exclusionDate(
    const DateTime& tradeDate, // when imnt thinks stock trade should happen
    const CAsset*   asset) {
    static const string method("PhysicalDelivery::exclusionDate");
    try {
        // what's going on here?
        // an imnt knows what the physical delivery is on T
        // we want the stock (or whatever) to drop out of the pricing
        // as of T+1 BEX (so EOD includes it, overnight does not)
        HolidayConstSP hols(AssetUtil::getHoliday(asset));
        DateTime tPlus1 = hols->addBusinessDays(tradeDate, 1);
        DateTime dropDate(tPlus1.getDate(), DateTime::BEFORE_EX_DIV_TIME);
        return dropDate;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** store results for the PHYSICAL_DELIVERY output request */
void PhysicalDelivery::recordPhysicalDelivery(
        Control*                     control,
        Results*                     results,    
        const string&                assetName,
        const PhysicalDeliveryArray* delivery) {
    static const string method = "PhysicalDelivery::recordPhysicalDelivery";
    try {
        OutputRequest* request = 
            control->requestsOutput(OutputRequest::PHYSICAL_DELIVERY);
        // if we don't add anything here then 'NotApplicable' will be added
        // by 'RiskManager' at the end of the pricing
        if (request && delivery && !delivery->empty()) {
            OutputNameSP name(new OutputName(assetName));
            IObjectSP    phys(copy(delivery));
            results->storeRequestResult(request, phys, name); 
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

void PhysicalDelivery::recordPhysicalDelivery(
    double          quantity, 
    double          price, 
    const DateTime& tradeDate, 
    const CAsset*   asset,
    Control*        control,
    Results*        results){
    static const string method = "PhysicalDelivery::recordPhysicalDelivery";
    try {
        OutputRequest* request = 
            control->requestsOutput(OutputRequest::PHYSICAL_DELIVERY);
        // if we don't add anything here then 'NotApplicable' will be added
        // by 'RiskManager' at the end of the pricing
        if (request) {
            const ICanPhysicallySettle* physical = dynamic_cast<const ICanPhysicallySettle*>(asset);
            if (!physical) {
                throw ModelException(method, 
                                     "asset (" + asset->getName() + 
                                     ") of type " + 
                                     asset->getClass()->getName() + 
                                     " can not physically settle");
            }

            PhysicalDeliveryByAssetArraySP byAsset(new PhysicalDeliveryByAssetArray(0));
            physical->delivery(tradeDate, quantity, price, byAsset.get());

            for (int i = 0; i < byAsset->size(); i++) {
                PhysicalDeliveryArray pda(1, *(*byAsset)[i].delivery.get());
                PhysicalDelivery::recordPhysicalDelivery(control,
                                                         results,
                                                         (*byAsset)[i].asset,
                                                         &pda);
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

PhysicalDelivery::PhysicalDelivery(double          quantity, 
                                   double          price, 
                                   const DateTime& tradeDate, 
                                   const DateTime& settleDate): 
    CObject(TYPE), quantity(quantity), price(price), 
    tradeDate(tradeDate), settleDate(settleDate) {}

PhysicalDelivery::PhysicalDelivery() : CObject(TYPE), quantity(0.0), price(0.0) {}

PhysicalDelivery::ByAsset::ByAsset(
    const string& asset, PhysicalDelivery* delivery): 
    CObject(TYPE), asset(asset), delivery(delivery) {}

PhysicalDelivery::ByAsset::ByAsset() : CObject(TYPE) {}

class PhysicalDeliveryHelper{
public:
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        clazz->setDescription("Container for physical settlement");
        REGISTER(PhysicalDelivery, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultPhysicalDelivery);
        FIELD(quantity, "quantity");
        FIELD(price, "price");
        FIELD(tradeDate, "trade date");
        FIELD(settleDate, "settle date");
    }
    static IObject* defaultPhysicalDelivery(){
        return new PhysicalDelivery();
    }

    static void loadByAsset(CClassSP& clazz){
        REGISTER(PhysicalDelivery::ByAsset, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultByAsset);
        FIELD(asset, "asset");
        FIELD(delivery, "delivery");
    }
    static IObject* defaultByAsset(){
        return new PhysicalDelivery::ByAsset();
    }
};

CClassConstSP const PhysicalDelivery::TYPE = CClass::registerClassLoadMethod(
    "PhysicalDelivery", typeid(PhysicalDelivery), PhysicalDeliveryHelper::load);

DEFINE_TEMPLATE_TYPE(PhysicalDeliveryArray);

CClassConstSP const PhysicalDelivery::ByAsset::TYPE = CClass::registerClassLoadMethod(
    "PhysicalDelivery::ByAsset", typeid(PhysicalDelivery::ByAsset), PhysicalDeliveryHelper::loadByAsset);

DEFINE_TEMPLATE_TYPE_WITH_NAME("PhysicalDelivery::ByAssetArray", PhysicalDeliveryByAssetArray);

static void loadICanPhysicallySettle(CClassSP& clazz){
    REGISTER_INTERFACE(ICanPhysicallySettle, clazz);
    EXTENDS(IObject);
}

CClassConstSP const ICanPhysicallySettle::TYPE = 
CClass::registerInterfaceLoadMethod("ICanPhysicallySettle", 
                                    typeid(ICanPhysicallySettle), 
                                    loadICanPhysicallySettle);

/** specialisations of arrayObjectCast */
/** Casts array element to an IObject */
IObjectConstSP arrayObjectCast<PhysicalDelivery>::toIObject(
    const PhysicalDelivery& value){
    IObjectConstSP objValue(IObjectConstSP::attachToRef(&value));
    return objValue;
}

/** Casts array element to an IObject */
IObjectSP arrayObjectCast<PhysicalDelivery>::toIObject(PhysicalDelivery& value){
    return PhysicalDeliverySP::attachToRef(&value);
}

/** Turns the IObjectSP into a PhysicalDelivery */
PhysicalDelivery arrayObjectCast<PhysicalDelivery>::fromIObject(IObjectSP& value){
    PhysicalDelivery *dtPtr = dynamic_cast<PhysicalDelivery *>(value.get());
    if (!dtPtr){
        throw ModelException("arrayObjectCast::fromIObject", "Object is not"
                             " a PhysicalDelivery");
    }
    return *dtPtr;
}

// explicit clone for arrays of PhysicalDeliverys - for performance
IObject* arrayClone<PhysicalDelivery>::clone(const CArray* arrayToClone){
    const PhysicalDeliveryArray& theArray = 
        static_cast<const PhysicalDeliveryArray&>(*arrayToClone);
    return new PhysicalDeliveryArray(theArray);
}

// now for the ByAsset flavours
/** Casts array element to an IObject */
IObjectConstSP arrayObjectCast<PhysicalDelivery::ByAsset>::toIObject(
    const PhysicalDelivery::ByAsset& value){
    IObjectConstSP objValue(IObjectConstSP::attachToRef(&value));
    return objValue;
}

/** Casts array element to an IObject */
IObjectSP arrayObjectCast<PhysicalDelivery::ByAsset>::toIObject(PhysicalDelivery::ByAsset& value){
    return PhysicalDeliveryByAssetSP::attachToRef(&value);
}

/** Turns the IObjectSP into a PhysicalDelivery::ByAsset */
PhysicalDelivery::ByAsset arrayObjectCast<PhysicalDelivery::ByAsset>::fromIObject(IObjectSP& value){
    PhysicalDelivery::ByAsset *dtPtr = dynamic_cast<PhysicalDelivery::ByAsset *>(value.get());
    if (!dtPtr){
        throw ModelException("arrayObjectCast::fromIObject", "Object is not"
                             " a PhysicalDelivery::ByAsset");
    }
    return *dtPtr;
}

// explicit clone for arrays of PhysicalDelivery::ByAssets - for performance
IObject* arrayClone<PhysicalDelivery::ByAsset>::clone(const CArray* arrayToClone){
    const PhysicalDeliveryByAssetArray& theArray = 
        static_cast<const PhysicalDeliveryByAssetArray&>(*arrayToClone);
    return new PhysicalDeliveryByAssetArray(theArray);
}


DRLIB_END_NAMESPACE

