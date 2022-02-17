//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CashFlow.hpp
//
//   Description : Well, it's a cashflow
//
//   Author      : Andrew J Swain
//
//   Date        : 7 February 2001
//
//----------------------------------------------------------------------------

#ifndef EDR_CASHFLOW_HPP
#define EDR_CASHFLOW_HPP

#include "edginc/DateTime.hpp"
#include "edginc/Expiry.hpp"

using namespace std;    // string

DRLIB_BEGIN_NAMESPACE
class CashFlow; // predeclaration in order to declare array
/** CashFlowArray is an array of CashFlow structures - not an array of
    pointers. Note that this approach forces the no argument CashFlow
    constructor to be public */
typedef array<CashFlow, CashFlow> CashFlowArray;
typedef smartPtr<CashFlowArray> CashFlowArraySP;
typedef smartConstPtr<CashFlowArray> CashFlowArrayConstSP;
#ifndef QLIB_CASHFLOW_CPP
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<CashFlow>);
EXTERN_TEMPLATE(class TOOLKIT_DLL array<CashFlow _COMMA_ CashFlow>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartConstPtr<CashFlowArray>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<CashFlowArray>);
EXTERN_TEMPLATE(IObjectSP TOOLKIT_DLL FieldGetInLine<CashFlowArray>(CashFlowArray* t));
EXTERN_TEMPLATE(void TOOLKIT_DLL FieldSetInLine<CashFlowArray>(CashFlowArray* t,
                                                   IObjectSP o));
#else
//// seem to need to instantiate certain templates before we get to cpp file
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<CashFlow>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL array<CashFlow _COMMA_ CashFlow>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartConstPtr<CashFlowArray>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<CashFlowArray>);
INSTANTIATE_TEMPLATE(IObjectSP TOOLKIT_DLL FieldGetInLine<CashFlowArray>(CashFlowArray* t));
INSTANTIATE_TEMPLATE(void TOOLKIT_DLL FieldSetInLine<CashFlowArray>(CashFlowArray* t,
                                                        IObjectSP o));
#endif

/** Defining a CashFlowCluster to be a CashFlowArrayArray. Note an array of
    structure not an array of pointers */
typedef array<CashFlowArray> CashFlowCluster;

class TOOLKIT_DLL CashFlow: public CObject {
public:
    friend class CashFlowHelper;
    static CClassConstSP const TYPE;

    CashFlow(const DateTime& date, double amount);
    CashFlow();
    ~CashFlow();

    DateTime date;
    double   amount;

    //// Returns true if dates are equal and if amounts are equal
    //// using Maths::isZero(. - .). Method added to support
    //// instantiating array template
    bool operator==(const CashFlow& rhs) const;

    /** Aggragates cashflows that occur at the same date and time. Modifies
        existing cashflow (so could end up with shorter array) */
    static void aggregate(CashFlowArray& cfArray);

    /** helper functions for sorting */
    static bool lessThenForDates(const CashFlow& x, const CashFlow& y);
    static bool lessThenForAmounts(const CashFlow& x, const CashFlow& y);

    /** Validates that dates within the cashflow are in increasing
        order. If failIfEmpty is true, then throws exception if cfArray
        is empty. The description string is used in the exception
        messages should describe what the cashflow represents */
    static void ensureDatesIncreasing(const CashFlowArray& cfArray,
                                      const string&        description,
                                      bool                 failIfEmpty);

    static CashFlowArraySP merge(const CashFlowArrayConstSP x, 
                                 const CashFlowArrayConstSP y);

    /** Simple Linear Interpolation routine on CashFlowArrays */
    static double interpolate(const CashFlowArray& cf,
                              const DateTime& date,
                              const bool extendFlat);

    /** This one requires an exact match on 'date' */
    static double interpolate(const CashFlowArray& cf,
                              const DateTime& date);

    /** Simple Linear Interpolation routine */
    static double interpolate(const ExpiryArray& exp,
                              const DoubleArray& vals,
                              const DateTime& today,
                              const DateTime& date,
                              const bool extendFlat);

    /** Returns the dates inside the CashFlowArray */
    static DateTimeArray dates(const CashFlowArray& cf);
    
    /** Returns the amounts inside the CashFlowArray */
    static CDoubleArraySP amounts(const CashFlowArray& cf);
    
    /** Constructs a CashFlowArray from dates and amounts arrays */
    static CashFlowArraySP createCashFlowArray(
        const DateTimeArray& dates,
        const CDoubleArray& amounts);
    
    /** Cooperate with strings */
    static CashFlowArray fromStringArray(const StringArray& asStrings);
    static StringArray toStringArray(const CashFlowArray& cfa);

    /** Accumulates the cashflows in the CashFlowArraySP passed in, so that 
        the cashflows in the returned CashFlowArraySP are:
           cf_out[i].amount = SUM_j{cf_in[j].amount} for all j <= i */
    static CashFlowArraySP accumulateCashFlows(CashFlowArraySP x);

private:
    class CFAddin;
    class GetCFAddin;
};

typedef smartPtr<CashFlow> CashFlowSP;

/** specialisations of arrayObjectCast */
template <> class TOOLKIT_DLL arrayObjectCast<CashFlow>{
public:
    /** Casts array element to an IObject */
    static IObjectConstSP toIObject(const CashFlow& value);

    /** Casts array element to an IObject */
    static IObjectSP toIObject(CashFlow& value);

    /** Turns the IObjectSP into a CashFlow */
    static const CashFlow& fromIObject(IObjectSP& value);
};

/** specialisation of arrayClone */
template <> class TOOLKIT_DLL arrayClone<CashFlow>{
public:
    /** Overridden for performance */
    static IObject* clone(const CArray* arrayToClone);
};

/** Defining a CashFlowCluster to be a CashFlowArrayArray. Note an array of
    structure not an array of pointers */
typedef array<CashFlowArray> CashFlowArrayArray;
// CashFlowArrayArray is the 'correct' name, cluster is just a convenience
typedef CashFlowArrayArray  CashFlowCluster;

/** specialisations of arrayObjectCast (needed as the array is not an array
    of pointers) */
template <> class TOOLKIT_DLL arrayObjectCast<CashFlowArray>{
public:
    /** Casts array element to an IObject */
    static IObjectConstSP toIObject(const CashFlowArray& value);

    /** Casts array element to an IObject */
    static IObjectSP toIObject(CashFlowArray& value);

    /** Turns the IObjectSP into a DateTime */
    static const CashFlowArray& fromIObject(IObjectSP& value);
};

/** specialisation of arrayClone */
template <> class TOOLKIT_DLL arrayClone<CashFlowArray>{
public:
    /** Overridden for performance */
    static IObject* clone(const CArray* arrayToClone);
};

typedef smartPtr<CashFlowCluster> CashFlowClusterSP;
typedef smartConstPtr<CashFlowCluster> CashFlowClusterConstSP;

/** wrapper around CashFlowArray that supports smart adding & scaling */
class TOOLKIT_DLL CashFlowList: public CObject,
                    virtual public CombinableResult {
public:
    static CClassConstSP const TYPE;

    CashFlowList(const CashFlowArray* cfl);

    /** scale by factor x (Implementation of CombinableResult) */
    virtual void scale(double x);
    /** add CashFlowList to this result (Implementation of
        CombinableResult) */
    virtual void add(const CombinableResult& x, double scaleFactor);

private:
    friend class CashFlowListHelper;
    CashFlowList();

    CashFlowArraySP cfl;
};

typedef smartPtr<CashFlowList> CashFlowListSP;

#ifndef QLIB_CASHFLOW_CPP
EXTERN_TEMPLATE(class TOOLKIT_DLL array<CashFlowArray _COMMA_ CashFlowArray>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<CashFlowList>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartConstPtr<CashFlowCluster>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<CashFlowCluster>);
#endif

DRLIB_END_NAMESPACE

#endif
