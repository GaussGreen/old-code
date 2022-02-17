//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : WrapperNameCollector.hpp
//
//   Description : Collects all wrapper names 
//
//   Author      : André Segger
//
//   Date        : 13 Jul 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDR_WRAPPERNAME_COLLECT_HPP
#define EDR_WRAPPERNAME_COLLECT_HPP
#include "edginc/Collector.hpp"
#include "edginc/AtomicArray.hpp"

DRLIB_BEGIN_NAMESPACE

/** This class should be scrapped. It is only used by one addin function
    (GET_WRAPPER_NAMES) which can probably be retired. Also the implementation
    is just daft - it requires every wrapper object to support it when that
    information could be read directly from the wrapper base class */
class RISKMGR_DLL WrapperNameCollector: public CObject, 
                            public virtual ICollector{

public:
    friend class WrapperNameCollectorHelper;
    static CClassConstSP const TYPE;

    /** default constructor */
    WrapperNameCollector();

    /** Add a name to the list */
    void addName(const string& assetName);

    /** triggers collection of all wrapper names */
    static CStringArraySP getWrapperNames(IObjectSP obj);

private:
    static void load(CClassSP& clazz);
    WrapperNameCollector(const WrapperNameCollector &rhs);
    WrapperNameCollector& operator=(const WrapperNameCollector& rhs);

    CStringArraySP     wrapperNames; // $unregistered
};

typedef smartPtr<WrapperNameCollector> WrapperNameCollectorSP;
typedef smartPtr<const WrapperNameCollector> WrapperNameCollectorConstSP;

DRLIB_END_NAMESPACE

#endif
