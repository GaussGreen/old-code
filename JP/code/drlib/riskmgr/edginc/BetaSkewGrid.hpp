//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : BetaSkewGrid.hpp
//
//   Description : Class to be able to aggregate BetaSkewGridResultArrays.
//                 Contains the actual result data array plus a method to 
//                 aggregate other objects of the same type
//
//   Author      : Jose Hilera
//
//   Date        : July 2005
//
//----------------------------------------------------------------------------

#ifndef QR_COMBINABLE_BETA_SKEW_RESULT_HPP
#define QR_COMBINABLE_BETA_SKEW_RESULT_HPP

#include "edginc/smartPtr.hpp"
#include "edginc/BetaSkewGridResult.hpp"
#include "edginc/CombinableResult.hpp"


DRLIB_BEGIN_NAMESPACE

class RISKMGR_DLL BetaSkewGrid :
    public CObject,
    virtual public CombinableResult
{
public:
    static CClassConstSP const TYPE;

    BetaSkewGrid(int size = 0);

    virtual ~BetaSkewGrid();

    void setPoint(int index, BetaSkewGridResultSP bs);

    /** scale by factor x (implementation of CombinableResult) */
    virtual void scale(double x);
    
    /** add another BetaSkewGrid object to this result
     * (implementation of CombinableResult) */
    virtual void add(const CombinableResult& x, double scaleFactor);

    /** write object out in 'output' format - ie suitable for comparing
     * regression files with */
    virtual void outputWrite(const string& linePrefix,
                             const string& prefix, 
                             ostream& stream) const;

private:
    // Reflection
    static void load(CClassSP& clazz);
    static IObject* defaultBetaSkewGrid();

    //Fields
    BetaSkewGridResultArraySP gridResult;  // The actual data results
};

typedef smartPtr<BetaSkewGrid> BetaSkewGridSP;
typedef smartConstPtr<BetaSkewGrid> BetaSkewGridConstSP;

DRLIB_END_NAMESPACE

#endif
