//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : AssetCcyCollector.hpp
//
//   Description : asset currency collector class
//
//   Author      : André Segger
//
//   Date        : 13 Jul 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDR_ASSETCCY_COLLECT_HPP
#define EDR_ASSETCCY_COLLECT_HPP
#include "edginc/Collector.hpp"

DRLIB_BEGIN_NAMESPACE
class YieldCurve;

/** A class to cross validate all fx rates used */
class MARKET_DLL AssetCcyCollector: public CObject,
                         public virtual ICollector {
public:
    friend class AssetCcyCollHelper;
    static CClassConstSP const TYPE;

    /** For use by 'accept' methods. Check's the currency is the same as that
        required. */ 
    void currencyValidate(const string& isoCode,
                          const string& assetName);

    /** Can be used to pass to 'accept' directly */
    AssetCcyCollector(const string& imntIsoCode);

    /** starts off empty, takes first ccy it finds as the one
        to validate subsequent ccy's against */
    AssetCcyCollector();

    /** triggers fx rate validation for a given object */
    static void validateAllCurrencies(IObjectConstSP    obj, 
                                      const YieldCurve* discCcy);

    /** triggers fx rate validation for a given object */
    static void validateAllCurrencies(IObjectConstSP obj,
                                      const string& discCcy);

    /** what ccy are we validating for ? */
    const string& ccyCode() const;

private:
    static void load(CClassSP& clazz);
    AssetCcyCollector(const AssetCcyCollector &rhs);
    AssetCcyCollector& operator=(const AssetCcyCollector& rhs);
    string           imntIsoCode; // $unregistered
    bool             empty; // $unregistered
};

typedef smartPtr<AssetCcyCollector> AssetCcyCollectorSP;
typedef smartPtr<const AssetCcyCollector> AssetCcyCollectorConstSP;

DRLIB_END_NAMESPACE

#endif
