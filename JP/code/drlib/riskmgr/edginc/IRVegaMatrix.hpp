//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives
//
//   Filename    : IRVegaMatrix.hpp
//
//   Description : IR vega pointwise sensitivity
//                 This is essentially the same of IRVegaPointwise but
//                 changed to store the array of results inside an object
//                 rather than storing the actual array.
//                 IRVegaPointwise is now deprecated.
//
//   Author      : Jose Hilera
//
//   Cautions    : The sensitivity implemented here is derived from 
//                 IRVegaPointwise, which is to be removed. To avoid duplication
//                 IRVegaPointwise::IShift and IRVegaPointwise::IRestorableShift
//                 are still used (and, indirectly, also 
//                 IRVegaPointwise::ISensitivePoints) - they should be moved into
//                 this class when IRVegaPointwise is finally removed, along the 
//                 methods in IRVegaPointwise which this class is currently 
//                 inheriting (NOT calculate, which VegaProxyMatrix overrides!)
//
//   Date        : July 2005
//
//----------------------------------------------------------------------------

#ifndef IRVEGAMATRIX_HPP
#define IRVEGAMATRIX_HPP

#include "edginc/IRGridShift.hpp"
#include "edginc/Additive.hpp"
#include "edginc/IRVegaPointwise.hpp"


DRLIB_BEGIN_NAMESPACE


class VegaProxyMatrix;

/** Sens Control for matrix IR vega */
class RISKMGR_DLL IRVegaMatrix: public IRVegaPointwise {
public:
    static CClassConstSP const TYPE;
    const static string NAME;
    const static double DEFAULT_SHIFT;

    /** constructor with explicit shift size */
    IRVegaMatrix(double shiftSize);

    /** implements a one sided vector derivative for each instance of the
        market data which is sensitive to this SensControl */
    virtual void calculate(TweakGroup*  tweakGroup,
                           CResults*    results);

private:
    friend class IRVegaMatrixHelper;
    /** for reflection */
    IRVegaMatrix();
    IRVegaMatrix(const IRVegaMatrix &rhs);
    IRVegaMatrix& operator=(const IRVegaMatrix& rhs);
};

typedef smartConstPtr<IRVegaMatrix> IRVegaMatrixConstSP;
typedef smartPtr<IRVegaMatrix> IRVegaMatrixSP;


DRLIB_END_NAMESPACE

#endif
