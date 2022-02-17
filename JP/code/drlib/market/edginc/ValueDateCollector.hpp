//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ValueDateCollector.hpp
//
//   Description : Value date collector class
//
//   Author      : Andre Segger
//
//   Date        : 30 Apr 2001
//
//
//----------------------------------------------------------------------------

#ifndef CVALUEDATE_COLLECT_HPP
#define CVALUEDATE_COLLECT_HPP
#include "edginc/Collector.hpp"
#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE

/** A class to validate cross validate start dates */
class MARKET_DLL CValueDateCollector: public CObject, 
                           public virtual ICollector {
public:
    friend class CValueDateCollHelper;
    static CClassConstSP const TYPE;
    void valueDateValidate(const DateTime& assetValueDate,
                           const string&   text);

    CValueDateCollector(const DateTime& valueDateToCompare,
                        const string& source);

    static void validateAll(IObjectSP       obj,
                            const DateTime& valueDate,
                            const string&   text);

private:
    static void load(CClassSP& clazz);
    CValueDateCollector(const CValueDateCollector &rhs);
    CValueDateCollector& operator=(const CValueDateCollector& rhs);

    DateTime            valueDate; // $unregistered
    string              source; // $unregistered
};

typedef smartPtr<CValueDateCollector> CValueDateCollectorSP;
typedef smartPtr<const CValueDateCollector> CValueDateCollectorConstSP;

DRLIB_END_NAMESPACE

#endif
