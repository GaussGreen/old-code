//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FutureExpiryCollector.hpp
//
//   Description : Future expiry date collector class
//
//   Author      : Andre Segger
//
//   Date        : 30 Apr 2001
//
//
//----------------------------------------------------------------------------

#ifndef FUTURE_COLLECT_HPP
#define FUTURE_COLLECT_HPP
#include <string>
#include "edginc/smartPtr.hpp"
#include "edginc/Collector.hpp"
#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE

/** A class to validate cross validate start dates */
class MARKET_DLL FutureExpiryCollector: public CObject, public virtual ICollector {
public:
    static CClassConstSP const TYPE;
    void expiryValidate(const DateTime& futureExpiryDate,
                        const string&   source);
    FutureExpiryCollector(const  DateTime& optionMaturity);
private:
    static void load(CClassSP& clazz);
    FutureExpiryCollector(const FutureExpiryCollector &rhs);
    FutureExpiryCollector& operator=(const FutureExpiryCollector& rhs);

    DateTime    maturityDate; // $unregistered
};

typedef smartPtr<FutureExpiryCollector> FutureExpiryCollectorSP;
typedef smartPtr<const FutureExpiryCollector> FutureExpiryCollectorConstSP;

DRLIB_END_NAMESPACE

#endif
