//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PhysicalDelivery.hpp
//
//   Description : Results data holder for what needs to be physically settled
//
//   Author      : Andrew J Swain
//
//   Date        : 26 February 2004
//
//
//----------------------------------------------------------------------------

#ifndef PHYSICALDELIVERY_HPP
#define PHYSICALDELIVERY_HPP

#include "edginc/DateTime.hpp"
#include "edginc/Asset.hpp"
#include "edginc/Control.hpp"
#include "edginc/Results.hpp"

DRLIB_BEGIN_NAMESPACE

class PhysicalDelivery;
typedef array<PhysicalDelivery>              PhysicalDeliveryArray;
typedef smartPtr<PhysicalDeliveryArray>      PhysicalDeliveryArraySP;
typedef smartConstPtr<PhysicalDeliveryArray> PhysicalDeliveryArrayConstSP;
typedef smartPtr<PhysicalDelivery>           PhysicalDeliverySP;
typedef smartConstPtr<PhysicalDelivery>      PhysicalDeliveryConstSP;

class MARKET_DLL PhysicalDelivery: public CObject {
public:
    static CClassConstSP const TYPE;

    // when physical should drop out of instrument price
    static DateTime exclusionDate(
        const DateTime& tradeDate, // when imnt thinks stock trade should happen
        const CAsset*   asset);

    /** store results for the PHYSICAL_DELIVERY output request */
    static void recordPhysicalDelivery(
        Control*                     control,
        Results*                     results,    
        const string&                assetName,
        const PhysicalDeliveryArray* delivery);

    static void recordPhysicalDelivery(
        double          quantity, 
        double          price, 
        const DateTime& tradeDate, 
        const CAsset*   asset,
        Control*        control,
        Results*        results);

    PhysicalDelivery(double          quantity, 
                     double          price, 
                     const DateTime& tradeDate, 
                     const DateTime& settleDate);

    /** Public default constructor to allow creation of PhysicalDeliveryArray's */ 
    PhysicalDelivery();

    // inner class - for use by ICanPhysicallySettle
    class MARKET_DLL ByAsset: public CObject {
    public:
        static CClassConstSP const TYPE;
        string             asset;
        PhysicalDeliverySP delivery;

        ByAsset(const string& asset, PhysicalDelivery* delivery);

        ByAsset();
    };


private:
    friend class PhysicalDeliveryHelper;

    // fields
    double   quantity;    // how much to deliver
    double   price;       // at this price
    DateTime tradeDate;   // traded here
    DateTime settleDate;  // and delivered here
};

typedef smartPtr<PhysicalDelivery::ByAsset>      PhysicalDeliveryByAssetSP;
typedef smartConstPtr<PhysicalDelivery::ByAsset> PhysicalDeliveryByAssetConstSP;
typedef array<PhysicalDelivery::ByAsset>         PhysicalDeliveryByAssetArray;
typedef smartPtr<PhysicalDeliveryByAssetArray>   PhysicalDeliveryByAssetArraySP;

// interface for assets that can physically deliver themselves
class MARKET_DLL ICanPhysicallySettle: public virtual IObject {
public:
    static CClassConstSP const TYPE;

    // for given trade date, return settle date and name of each asset
    virtual void delivery(const DateTime&               tradeDate, 
                          double                        quantity, 
                          double                        price,
                          PhysicalDeliveryByAssetArray* byAsset) const = 0;
};

/** specialisations of arrayObjectCast */
template <> class MARKET_DLL arrayObjectCast<PhysicalDelivery>{
public:
    /** Casts array element to an IObject */
    static IObjectConstSP toIObject(const PhysicalDelivery& value);

    /** Casts array element to an IObject */
    static IObjectSP toIObject(PhysicalDelivery& value);

    /** Turns the IObjectSP into a PhysicalDelivery */
    static PhysicalDelivery fromIObject(IObjectSP& value);
};

/** specialisation of arrayClone */
template <> class MARKET_DLL arrayClone<PhysicalDelivery>{
public:
    /** Overridden for performance */
    static IObject* clone(const CArray* arrayToClone);
};

/** specialisations of arrayObjectCast */
template <> class MARKET_DLL arrayObjectCast<PhysicalDelivery::ByAsset>{
public:
    /** Casts array element to an IObject */
    static IObjectConstSP toIObject(const PhysicalDelivery::ByAsset& value);

    /** Casts array element to an IObject */
    static IObjectSP toIObject(PhysicalDelivery::ByAsset& value);

    /** Turns the IObjectSP into a PhysicalDelivery::ByAsset */
    static PhysicalDelivery::ByAsset fromIObject(IObjectSP& value);
};

/** specialisation of arrayClone */
template <> class MARKET_DLL arrayClone<PhysicalDelivery::ByAsset>{
public:
    /** Overridden for performance */
    static IObject* clone(const CArray* arrayToClone);
};

DRLIB_END_NAMESPACE

#endif
