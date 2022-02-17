/**
 * @file CorrSwapBasisAdjParallel.hpp
 */

#ifndef QLIB_CorrSwapBasisAdjParallel_H
#define QLIB_CorrSwapBasisAdjParallel_H

#include "edginc/Additive.hpp"
#include "edginc/AllNamesRiskPropertySensitivity.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(CorrSwapBasisAdjParallel)

/**
 * A greek calculated by tweaking all the market names' corr swap basis
 * adjustments
 */

class RISKMGR_DLL CorrSwapBasisAdjParallel: public AllNamesRiskPropertySensitivity ,
                                public virtual Additive {

    static void load(CClassSP& clazz);
    Deriv deriv() const;

public:

    static CClassConstSP const TYPE;

    static const double DEFAULT_SHIFT;
    static const double SENSITIVITY_UNIT;
    static const string NAME;

    CorrSwapBasisAdjParallel(double shiftSize = DEFAULT_SHIFT);
    ~CorrSwapBasisAdjParallel();
};

typedef smartPtr<CorrSwapBasisAdjParallel> CorrSwapBasisAdjParallelSP;

DRLIB_END_NAMESPACE

#endif
