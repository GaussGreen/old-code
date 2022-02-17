//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : BetaSkewGridResult.hpp
//
//   Description : Result for pointwise type beta skew sensitivities
//
//   Author      : Antoine Gregoire
//
//   Date        : 17-Jun-2005
//
//
//   $Log: $
//
//----------------------------------------------------------------------------

#ifndef QLIB_BETA_SKEW_GRID_RESULT_HPP
#define QLIB_BETA_SKEW_GRID_RESULT_HPP

#include "edginc/BetaSkewGridPoint.hpp"
#include "edginc/CombinableResult.hpp"

DRLIB_BEGIN_NAMESPACE

/** Result for pointwise type beta skew sensitivities */
class RISKMGR_DLL BetaSkewGridResult: public CObject,
                          public virtual CombinableResult {
public:
    static CClassConstSP const TYPE;

    virtual ~BetaSkewGridResult();

    BetaSkewGridResult(const BetaSkewGridPoint& point, double result);

    /** returns the grid point associated with this result */
    BetaSkewGridPoint getGridPoint() const;

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
    
    /** add BetaSkewGridResult object to this result (Implementation of
        CombinableResult) */
    virtual void add(const CombinableResult& x, double scaleFactor);

private:
    // Reflection
    BetaSkewGridResult();
    static void load(CClassSP& clazz);
    static IObject* defaultBetaSkewGridResult();
    
    // Fields
    BetaSkewGridPoint point;
    double result;
};

typedef smartConstPtr<BetaSkewGridResult> BetaSkewGridResultConstSP;
typedef smartPtr<BetaSkewGridResult> BetaSkewGridResultSP;

typedef array<BetaSkewGridResultSP, BetaSkewGridResult> BetaSkewGridResultArray;
typedef smartPtr<BetaSkewGridResultArray> BetaSkewGridResultArraySP;
typedef smartConstPtr<BetaSkewGridResultArray> BetaSkewGridResultArrayConstSP;

DRLIB_END_NAMESPACE

#endif //QLIB_BETA_SKEW_GRID_RESULT_HPP
