//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ExpiryResult.hpp
//
//   Description : Captures result of a result at a given expiry, ie stores
//                 an expiry and a double
//
//   Author      : Mark A Robson
//
//   Date        : 5 March 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_EXPIRY_RESULT_HPP
#define EDG_EXPIRY_RESULT_HPP

#include "edginc/Expiry.hpp"
#include "edginc/CombinableResult.hpp"

using namespace std;    // string

DRLIB_BEGIN_NAMESPACE

// Establish array type for ExpiryResult::createExpiryResultArray()

class ExpiryResult;

/** ExpiryResultArray is an array of ExpiryResult structures - not an array of
    pointers. Note that this approach forces the no argument ExpiryResult
    constructor to be public */

typedef array<ExpiryResult> ExpiryResultArray;
typedef smartPtr<ExpiryResultArray> ExpiryResultArraySP;
typedef smartConstPtr<ExpiryResultArray> ExpiryResultArrayConstSP;

#ifndef QLIB_EXPIRYRESULT_CPP
EXTERN_TEMPLATE(class TOOLKIT_DLL array<ExpiryResult _COMMA_ ExpiryResult>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartConstPtr<ExpiryResultArray>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<ExpiryResultArray>);
#else
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL array<ExpiryResult _COMMA_ ExpiryResult>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartConstPtr<ExpiryResultArray>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<ExpiryResultArray>);
#endif

/** Captures result of a result at a given expiry, ie stores an expiry
    and a double */

class TOOLKIT_DLL ExpiryResult: public CObject,
                    public CombinableResult {
public:
    static CClassConstSP const TYPE;
    friend class ExpiryResultHelper;

    virtual ~ExpiryResult();

    /* for simplicity allow array of ExpiryResults to be an array of
       structures rather than an array of pointers - dictates public
       default constructor */
    ExpiryResult();

    ExpiryResult(const ExpiryConstSP& expiry, double result);

    /** returns the expiry associated with this result */
    ExpiryConstSP getExpiry() const;

    /** returns the double detailing this result */
    double getResult() const;

    /** adds the supplied value to this object */
    void addToResult(double value);

    //// Returns true if expiries are equal and if "result"s are equal
    //// using Maths::isZero(. - .). Method added to support
    //// instantiating array template
    bool operator==(const ExpiryResult& rhs) const;

    /** write object out in 'output' format - ie suitable for comparing
        regression files with */
    virtual void outputWrite(const string& linePrefix,
                             const string& prefix, ostream& stream) const;

    /** scale by factor x (Implementation of CombinableResult) */
    virtual void scale(double x);
    /** add ExpiryResult object to this result (Implementation of
        CombinableResult) */
    virtual void add(const CombinableResult& x, double scaleFactor);

    /** Constructs an ExpiryResultArray from dates and amounts arrays */
    static ExpiryResultArraySP createExpiryResultArray(
        const ExpiryArray& expiries,
        const CDoubleArray& amounts);

    /** Returns index of the expiry in supplied ExpiryResultArray.
        The supplied base date is used in the comparison so that it is 
        possible to search for eg MaturityPeriods in ExpiryResultArrays 
        were the expiries are BenchmarkDates */
    static int search(const DateTime& base,
                      const ExpiryConstSP expiry,
                      const ExpiryResultArraySP expiryResults);

    /** Reports if the supplied ExpiryArray is a subset of the expiries in the
        ExpiryResultArray - NB: not the other way around!
        The supplied base date is used in the comparison so that it is 
        possible to compare eg MaturityPeriods and ExpiryResultArrays 
        were the expiries are BenchmarkDates */
    static int isSubset(const DateTime& base,
                                      const ExpiryArrayConstSP part,
                                      const ExpiryResultArraySP full);

private:
    ExpirySP  expiry;
    double    result;
};



/** specialisations of arrayObjectCast - needed to support arrays of 
    ExpiryResult */
template <> class TOOLKIT_DLL arrayObjectCast<ExpiryResult>{
public:
    /** Casts array element to an IObject */
    static IObjectConstSP toIObject(const ExpiryResult& value);

    /** Casts array element to an IObject */
    static IObjectSP toIObject(ExpiryResult& value);

    /** Turns the IObjectSP into a DateTime */
    static const ExpiryResult& fromIObject(IObjectSP& value);
};

DRLIB_END_NAMESPACE
#endif
