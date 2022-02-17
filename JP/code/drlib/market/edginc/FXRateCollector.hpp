//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CFXRateCollector.hpp
//
//   Description : fx rate collector class
//
//   Author      : Andre Segger
//
//   Date        : 01 May 2001
//
//
//----------------------------------------------------------------------------

#ifndef CFXRATE_COLLECT_HPP
#define CFXRATE_COLLECT_HPP
#include "edginc/smartPtr.hpp"
#include "edginc/Collector.hpp"
#include <map>

DRLIB_BEGIN_NAMESPACE

/** A class to cross validate all fx rates used */
class MARKET_DLL CFXRateCollector: public CObject, 
                        public virtual ICollector {
public:
    friend class CFXRateCollHelper;
    static CClassConstSP const TYPE;

    ~CFXRateCollector();

    /** check an fx rate against preceding rates. Fails if any rates have 
        the same name yet different rates */ 
    void fxRateValidate(const string& fxName,
                        const double& fxRate,
                        const string& source);
    CFXRateCollector();

    /** triggers fx rate validation for a given object */
    static void validateAllFXRates(IObjectSP obj);

private:
    static void load(CClassSP& clazz);
    CFXRateCollector(const CFXRateCollector &rhs);
    CFXRateCollector& operator=(const CFXRateCollector& rhs);

    /** the fx rate map stores a set of fx rates (names) and corresponding
        fx rates */
    map<string, double> fxRateMap; // $unregistered
};

typedef smartPtr<CFXRateCollector> CFXRateCollectorSP;
typedef smartPtr<const CFXRateCollector> CFXRateCollectorConstSP;

DRLIB_END_NAMESPACE

#endif
