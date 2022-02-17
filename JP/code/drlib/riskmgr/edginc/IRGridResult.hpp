//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IRGridResult.hpp
//
//   Description : Result for pointwise type IR vega sensitivities
//
//   Author      : Andrew J Swain
//
//   Date        : 26 February 2002
//
//
//----------------------------------------------------------------------------

#ifndef IRGRIDRESULT_HPP
#define IRGRIDRESULT_HPP

#include "edginc/IRGridPointAbs.hpp"
#include "edginc/CombinableResult.hpp"

using namespace std;    // string

DRLIB_BEGIN_NAMESPACE

/** Result for pointwise type IR vega sensitivities */
class RISKMGR_DLL IRGridResult: public CObject,
                    public CombinableResult {
public:
    static CClassConstSP const TYPE;
    friend class IRGridResultHelper;

    virtual ~IRGridResult();

    /* for simplicity allow array of IRGridResults to be an array of
       structures rather than an array of pointers - dictates public
       default constructor */
    IRGridResult();

    IRGridResult(IRGridPointAbsConstSP point, double result);

    IRGridResult(IRGridPointSP point, double result);

    /** returns the grid point associated with this result */
    IRGridPointSP getGridPoint() const;

    /** returns the double detailing this result */
    double getResult() const;

    /** adds the supplied value to this object */
    void addToResult(double value);

    /** write object out in 'output' format - ie suitable for comparing
        regression files with */
    virtual void outputWrite(const string& linePrefix,
                             const string& prefix, ostream& stream) const;

    /** scale by factor x (Implementation of CombinableResult) */
    virtual void scale(double x);
    /** add IRGridResult object to this result (Implementation of
        CombinableResult) */
    virtual void add(const CombinableResult& x, double scaleFactor);

private:
    IRGridPointSP point;
    double        result;
};



/** specialisations of arrayObjectCast - needed to support arrays of 
    IRGridResult */
template <> class RISKMGR_DLL arrayObjectCast<IRGridResult>{
public:
    /** Casts array element to an IObject */
    static IObjectConstSP toIObject(const IRGridResult& value);

    /** Casts array element to an IObject */
    static IObjectSP toIObject(IRGridResult& value);

    /** Turns the IObjectSP into a DateTime */
    static IRGridResult fromIObject(IObjectSP& value);
};

/** IRGridResultArray is an array of IRGridResult structures - not an array of
    pointers. Note that this approach forces the no argument IRGridResult
    constructor to be public */
typedef array<IRGridResult, IRGridResult> IRGridResultArray;
typedef smartPtr<IRGridResultArray> IRGridResultArraySP;
typedef smartConstPtr<IRGridResultArray> IRGridResultArrayConstSP;

DRLIB_END_NAMESPACE
#endif
