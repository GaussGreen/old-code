/**
 * @file CreditIndexSpreadPointwise.hpp
 */

#ifndef QLIB_CreditIndexSpreadPointwise_H
#define QLIB_CreditIndexSpreadPointwise_H

#include "edginc/ExpiryWindow.hpp"
#include "edginc/RiskProperty.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Tag class denoting "index spread at an expiry" as a property of market names.
 *
 * This class doesn't have any methods or members --- it just acts as a "tag"
 * for instantiating a family of templates,
 *
 * See VolPointwise for a fully documented analogous example.
 */

struct CreditIndexSpreadPointwise: CObject {

    static CClassConstSP const TYPE;
    CreditIndexSpreadPointwise(); ~CreditIndexSpreadPointwise();

    /**
     * CreditIndexSpreadPointwise exists in multiple instances for each name, and requires an
     * ExpiryWindow for full specification.
     */
    typedef ExpiryWindow Qualifier;

    enum {
        /**
         * CreditIndexSpreadPointwise is continuous, i.e. tweaks to it can be made
         * arbitrarily small.
         *
         * This is carried through to Sensitivity::discreteShift()
         * via RiskProperty<PROPERTY>::discrete().
         */
        discrete = 0
    };
};

//declaration of template specialisation
template <>
RiskMappingMatrixConstSP RiskProperty<CreditIndexSpreadPointwise>::riskMappingMatrix(
    IObjectConstSP world,
    OutputNameConstSP name) const;

DRLIB_END_NAMESPACE
#endif
