//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : BarrierLevel.hpp
//
//   Description : Results data holder for barrier level today
//
//   Author      : Andrew J Swain
//
//   Date        : 10 June 2003
//
//
//----------------------------------------------------------------------------

#ifndef BARRIER_LEVEL__HPP
#define BARRIER_LEVEL__HPP

#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE

class UTIL_DLL BarrierLevel: public CObject {
public:
    static CClassConstSP const TYPE;

    // how many days forward we report barriers
    static const int DATE_WINDOW;

    // temporary defaulting until all products catch up with change
    BarrierLevel(bool            isUp, 
                 const DateTime& date, 
                 double          level, 
                 bool            isContinuous = true);

    // given today's date, when should we report barriers up to ?
    static DateTime barrierWindow(const DateTime& today);
    
    /** Public default constructor to allow creation of BarrierLevelArray's */ 
    BarrierLevel();

    //// Routes through equalTo. Method added to support
    //// instantiating array template
    bool operator==(const BarrierLevel& rhs) const;
private:
    friend class BarrierLevelHelper;

    // fields
    bool     isUp;
    DateTime date;
    double   level;
    bool     isContinuous;
};

/** specialisations of arrayObjectCast */
template <> class UTIL_DLL arrayObjectCast<BarrierLevel>{
public:
    /** Casts array element to an IObject */
    static IObjectConstSP toIObject(const BarrierLevel& value);

    /** Casts array element to an IObject */
    static IObjectSP toIObject(BarrierLevel& value);

    /** Turns the IObjectSP into a BarrierLevel */
    static BarrierLevel fromIObject(IObjectSP& value);
};

/** specialisation of arrayClone */
template <> class UTIL_DLL arrayClone<BarrierLevel>{
public:
    /** Overridden for performance */
    static IObject* clone(const CArray* arrayToClone);
};

typedef smartPtr<BarrierLevel> BarrierLevelSP;
typedef smartConstPtr<BarrierLevel> BarrierLevelConstSP;
#ifndef QLIB_BARRIERLEVEL_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartConstPtr<BarrierLevel>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<BarrierLevel>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartConstPtr<BarrierLevel>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<BarrierLevel>);
#endif

typedef array<BarrierLevel>              BarrierLevelArray;
typedef smartPtr<BarrierLevelArray>      BarrierLevelArraySP;
typedef smartConstPtr<BarrierLevelArray> BarrierLevelArrayConstSP;
#ifndef QLIB_BARRIERLEVEL_CPP
EXTERN_TEMPLATE(class UTIL_DLL array<BarrierLevel>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartConstPtr<BarrierLevelArray>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<BarrierLevelArray>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL array<BarrierLevel>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartConstPtr<BarrierLevelArray>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<BarrierLevelArray>);
#endif

DRLIB_END_NAMESPACE

#endif
