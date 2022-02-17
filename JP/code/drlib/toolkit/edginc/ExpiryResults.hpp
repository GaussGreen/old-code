//----------------------------------------------------------------------------
//
//   Group       : GCCT Derivatives Research
//
//   Filename    : ExpiryResults.hpp
//
//   Description : Captures double results of a result at a given expiry, ie stores
//                 an expiry and a double array. Strickly based on ExpiryResults
//
//   Author      : Sean Chen
//
//   Date        : 17 May 2005
//
//
//----------------------------------------------------------------------------

#ifndef EDG_EXPIRY_RESULTS_HPP
#define EDG_EXPIRY_RESULTS_HPP

#include "edginc/Expiry.hpp"
#include "edginc/CombinableResult.hpp"

using namespace std;    // string

DRLIB_BEGIN_NAMESPACE

/** Captures results of a result at a given expiry, ie stores an expiry
    and a double array */

class TOOLKIT_DLL ExpiryResults: public CObject,
                    public CombinableResult {
public:
    static CClassConstSP const TYPE;
    friend class ExpiryResultsHelper;

    virtual ~ExpiryResults();

    /* for simplicity allow array of ExpiryResults to be an array of
       structures rather than an array of pointers - dictates public
       default constructor */
    ExpiryResults();

    ExpiryResults(const ExpiryConstSP& expiry, const CDoubleArray& result);

    /** returns the expiry associated with this result */
    ExpiryConstSP getExpiry() const;

    /** returns the double detailing this result */
    CDoubleArray getResult() const;

    /** adds the supplied value to this object */
    void addToResult(const CDoubleArray& value);

    /** write object out in 'output' format - ie suitable for comparing
        regression files with */
    virtual void outputWrite(const string& linePrefix,
                             const string& prefix, ostream& stream) const;

    /** scale by factor x (Implementation of CombinableResult) */
    virtual void scale(double x);
    /** add ExpiryResults object to this result (Implementation of
        CombinableResult) */
    virtual void add(const CombinableResult& x, double scaleFactor);

private:
    ExpirySP  expiry;
    CDoubleArray    result;
};



/** specialisations of arrayObjectCast - needed to support arrays of 
    ExpiryResults */
template <> class TOOLKIT_DLL arrayObjectCast<ExpiryResults>{
public:
    /** Casts array element to an IObject */
    static IObjectConstSP toIObject(const ExpiryResults& value);

    /** Casts array element to an IObject */
    static IObjectSP toIObject(ExpiryResults& value);

    /** Turns the IObjectSP into a DateTime */
    static ExpiryResults fromIObject(IObjectSP& value);
};

/** ExpiryResultsArray is an array of ExpiryResults structures - not an array of
    pointers. Note that this approach forces the no argument ExpiryResults
    constructor to be public */
typedef array<ExpiryResults> ExpiryResultsArray;
typedef smartPtr<ExpiryResultsArray> ExpiryResultsArraySP;
typedef smartConstPtr<ExpiryResultsArray> ExpiryResultsArrayConstSP;

DRLIB_END_NAMESPACE
#endif
