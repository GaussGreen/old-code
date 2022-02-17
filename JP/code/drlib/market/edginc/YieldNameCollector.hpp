//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : YieldNameCollector.hpp
//
//   Description : Collects all yield names 
//
//   Author      : Stephen Hope
//
//   Date        : 23 Jan 2002
//
//
//----------------------------------------------------------------------------

#ifndef EDR_YIELDNAME_COLLECT_HPP
#define EDR_YIELDNAME_COLLECT_HPP
#include "edginc/Collector.hpp"
#include "edginc/AtomicArray.hpp"


DRLIB_BEGIN_NAMESPACE

/** A class to collect all yield curve names.
    This class should be scrapped. It is only used by one addin function
    (GET_YIELD_NAMES) which can probably be retired. Also the implementation
    is just daft - it requires each yield curve object to support it when that
    information could be read directly from the yield curve [base] class */
class MARKET_DLL YieldNameCollector: public CObject, 
                          public virtual ICollector {
public:
    friend class YieldNameCollectorHelper;
    static CClassConstSP const TYPE;

    /** default constructor */
    YieldNameCollector();

    /** Implementation of IAction interface */
    bool invoke(const IObjectSP& obj);

    /** Implementation of IAction interface */
    void addName(const string& assetName);

    /** triggers collection of all Yield names */
    static CStringArraySP getYieldNames(IObjectSP obj);

private:
    static void load(CClassSP& clazz);
    YieldNameCollector(const YieldNameCollector &rhs);
    YieldNameCollector& operator=(const YieldNameCollector& rhs);

    CStringArraySP     yieldNames; // $unregistered
};

typedef smartPtr<YieldNameCollector> YieldNameCollectorSP;
typedef smartPtr<const YieldNameCollector> YieldNameCollectorConstSP;

DRLIB_END_NAMESPACE

#endif
