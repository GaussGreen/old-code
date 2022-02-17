/**
 * @file NamedRiskQuantity.hpp
 */

#ifndef DRLIB_NamedRiskQuantity_H
#define DRLIB_NamedRiskQuantity_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/Results_forward.hpp"
#include "edginc/RiskQuantity.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(IResultsIdentifier)
FORWARD_DECLARE(NamedRiskQuantity)

/**
 * A RiskQuantity which knows where it is to be stored in the Results.
 *
 * Represents e.g. "please estimate IBM 3M vega by a 1bp one-sided tweak
 * and store the result under 'VEGA_POINTWISE', 'IBM', '3M' scaled into
 * 1bp units".
 *
 * See RiskQuantity, IResultsIdentifier.  It's convenient to have a class
 * bringing these together because there is a little bit of logic in
 * storeResult().
 *
 * Instances of this class are created in e.g.
 * PerNameRiskPropertySensitivity::nameRiskQuantities().
 *
 * For an overview of the "declarative" sensitivities framework of which these
 * classes are a part, see IRiskQuantityFactory.
 */

class RISKMGR_DLL NamedRiskQuantity: public CObject {

    static void load(CClassSP& clazz);
    NamedRiskQuantity(const NamedRiskQuantity& rhs);
    NamedRiskQuantity& operator=(const NamedRiskQuantity& rhs);

public:

    static CClassConstSP const TYPE;

    /**
     * Constructor.
     *
     * @param unit   Typically 1, or e.g. 0.0001 to get Results entry
     *               scaled in basis points rather than natural units
     */

    NamedRiskQuantity(RiskQuantityConstSP riskQuantity,
                      IResultsIdentifierConstSP resultsName,
                      double unit = 1.);

    /**
     * Constructor returning smartPtr.
     */

    static NamedRiskQuantitySP SP(RiskQuantityConstSP riskQuantity,
                                  IResultsIdentifierConstSP resultsName,
                                  double unit = 1.);

    ~NamedRiskQuantity();

    RiskQuantityConstSP const riskQuantity;
    IResultsIdentifierConstSP const resultsName;
    double unit;

    /**
     * Evaluate the RiskQuantity and store its value in the Results under the
     * IResultsIdentifier, scaled by unit.
     *
     * Calls riskQuantity.value(vals, dsts) to get the estimated greek, and
     * then resultsName.storeResult(results, value) to put it away.
     */

    virtual void storeResult(const CDoubleArray& vals, const CDoubleArray& dsts,
                             CResultsSP results) const;
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
