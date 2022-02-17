/**
 * @file ExpiryPairResult.hpp
 */

#ifndef QLIB_ExpiryPairResult_H
#define QLIB_ExpiryPairResult_H

#include "edginc/Expiry.hpp"
#include "edginc/ExpiryPair.hpp"
#include "edginc/CombinableResult.hpp"

using namespace std;    // string

DRLIB_BEGIN_NAMESPACE

// Establish array type for ExpiryPairResult::createExpiryPairResultArray()

class ExpiryPairResult;

/** ExpiryPairResultArray is an array of ExpiryPairResult structures - not an array of
    pointers. Note that this approach forces the no argument ExpiryPairResult
    constructor to be public */

typedef array<ExpiryPairResult> ExpiryPairResultArray;
typedef smartPtr<ExpiryPairResultArray> ExpiryPairResultArraySP;
typedef smartConstPtr<ExpiryPairResultArray> ExpiryPairResultArrayConstSP;

#ifndef QLIB_EXPIRYPAIRRESULT_CPP
EXTERN_TEMPLATE(class TOOLKIT_DLL array<ExpiryPairResult _COMMA_ ExpiryPairResult>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartConstPtr<ExpiryPairResultArray>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<ExpiryPairResultArray>);
#else
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL array<ExpiryPairResult _COMMA_ ExpiryPairResult>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartConstPtr<ExpiryPairResultArray>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<ExpiryPairResultArray>);
#endif

/** Captures result of a result at a given expiry, ie stores an expiry
    and a double */

class TOOLKIT_DLL ExpiryPairResult: public CObject,
                    public CombinableResult {
public:
    static CClassConstSP const TYPE;
    friend class ExpiryPairResultHelper;

    virtual ~ExpiryPairResult();

    /* for simplicity allow array of ExpiryPairResults to be an array of
       structures rather than an array of pointers - dictates public
       default constructor */
    ExpiryPairResult();

    ExpiryPairResult(const ExpiryPairConstSP& expiryPair, double result);

    /** returns the expiry associated with this result */
    ExpiryPairConstSP getExpiryPair() const;

    /** returns the double detailing this result */
    double getResult() const;

    /** adds the supplied value to this object */
    void addToResult(double value);

    //// Returns true if expiries are equal and if "result"s are equal
    //// using Maths::isZero(. - .). Method added to support
    //// instantiating array template
    bool operator==(const ExpiryPairResult& rhs) const;

    /** write object out in 'output' format - ie suitable for comparing
        regression files with */
    virtual void outputWrite(const string& linePrefix,
                             const string& prefix, ostream& stream) const;

    /** scale by factor x (Implementation of CombinableResult) */
    virtual void scale(double x);
    /** add ExpiryPairResult object to this result (Implementation of
        CombinableResult) */
    virtual void add(const CombinableResult& x, double scaleFactor);

    /** Constructs an ExpiryPairResultArray from dates and amounts arrays */
    static ExpiryPairResultArraySP createExpiryPairResultArray(
        const ExpiryPairArray& expiryPairs,
        const CDoubleArray& amounts);

    /** Returns index of the expiry in supplied ExpiryPairResultArray.
        The supplied base date is used in the comparison so that it is 
        possible to search for eg MaturityPeriods in ExpiryPairResultArrays 
        were the expiries are BenchmarkDates */
    static int search(const DateTime& base,
                      const ExpiryPairConstSP expiryPair,
                      const ExpiryPairResultArraySP expiryResults);

    /** Reports if the supplied ExpiryPairArray is a subset of the expiries in the
        ExpiryPairResultArray - NB: not the other way around!
        The supplied base date is used in the comparison so that it is 
        possible to compare eg MaturityPeriods and ExpiryPairResultArrays 
        were the expiries are BenchmarkDates */
    static int isSubset(const DateTime& base,
                                      const ExpiryPairArrayConstSP part,
                                      const ExpiryPairResultArraySP full);

private:
    ExpiryPairSP  expiryPair;
    double    result;
};



/** specialisations of arrayObjectCast - needed to support arrays of 
    ExpiryPairResult */
template <> class TOOLKIT_DLL arrayObjectCast<ExpiryPairResult>{
public:
    /** Casts array element to an IObject */
    static IObjectConstSP toIObject(const ExpiryPairResult& value);

    /** Casts array element to an IObject */
    static IObjectSP toIObject(ExpiryPairResult& value);

    /** Turns the IObjectSP into a DateTime */
    static const ExpiryPairResult& fromIObject(IObjectSP& value);
};

DRLIB_END_NAMESPACE
#endif
