/**
 * @file VegaParallel2Sided.cpp
 */

#ifndef QLIB_VegaParallel2Sided_H
#define QLIB_VegaParallel2Sided_H

#include "edginc/config.hpp"
#include "edginc/config.hpp"
#include "edginc/MultiTweakGroup.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/PropertyRiskAxis.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/PropertyTweakHypothesis.hpp"
#include "edginc/RiskQuantity.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/NamedRiskQuantity.hpp"
#include "edginc/GenericScalarTwoSidedShift.hpp"
#include "edginc/VolParallel.hpp"
#include "edginc/CreditTweak.hpp"

DRLIB_BEGIN_NAMESPACE
FORWARD_DECLARE(VegaParallel2Sided)
//// For windows dlls, the specific GenericScalarTwoSidedShift must be
//// externed/instantiated below in order to see the static data member in
//// another dll
#ifndef QLIB_VEGAPARALLEL2SIDED_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL GenericScalarTwoSidedShift<VolParallel _COMMA_
                false _COMMA_ true>);
#else
INSTANTIATE_TEMPLATE(class RISKMGR_DLL GenericScalarTwoSidedShift<VolParallel 
                     _COMMA_ false _COMMA_ true>);
#endif

/**
 * Calculates first and second derivatives of instrument price w.r.t. bulk
 * shifts in vol
 */

class RISKMGR_DLL VegaParallel2Sided:
    public GenericScalarTwoSidedShift<VolParallel, false, true> {

    typedef GenericScalarTwoSidedShift<VolParallel, false, true> Super;

    static void load(CClassSP& clazz);

    friend class DDeltaDVol;
    VegaParallel2Sided(double shiftSize, IModel* model, CControl* control);

    IObjectSP VegaPar; // field was always spurious---kept for back-compat

    static IObject* defaultConstructor();
public:
    VegaParallel2Sided(double shiftSize);

    static CClassConstSP const TYPE;

};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
