/**
 * @file RiskAxis.hpp
 */

#ifndef DRLIB_RiskAxis_H
#define DRLIB_RiskAxis_H

#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(RiskAxis)
FORWARD_DECLARE(IRiskAxis)

/**
 * "Frozen" representation of an IRiskAxis, for database storage.
 *
 * To implement RiskMapping we need to be able to store a large number of
 * different kinds of IRiskAxis in Pyramid.  Rather than hardcode them all into
 * the Pyramid data model, we've introduced a small number of factory classes
 * which store the internal QLib class name of an IRiskAxis together with its
 * constructor arguments, and use them to implement the thawed() method
 * defined here.
 *
 * This should really be called FrozenRiskAxis to avoid confusing us, because
 * it's not a base class for IRiskAxis implementations ... but on the other
 * hand that would confuse the Pyramid guys.
 *
 * The mechanism is as follows:
 *
 *    -  Every IRiskAxis must implement IRiskAxis::frozen(), returning
 *       a RiskAxis "storage representation"
 *
 *    -  In many cases it's sufficient to use frozenSimple(), having
 *       said which fields on the IRiskAxis class carry the market data
 *       name and qualifier (if any) at "load" time, using registerSimple().
 *       Example: PropertyRiskAxis<PROPERTY>::frozen().
 *
 *    -  In other cases it's easiest not to bother with a separate
 *       frozen representation: then the IRiskAxis can itself implement
 *       RiskAxis, with its IRiskAxis::frozen() and RiskAxis::thawed() both
 *       just returning "this". Example: CalibratorFieldRiskAxis::frozen().
 *
 * For an overview of the "declarative" sensitivities framework of which
 * IRiskAxis is a part, see IRiskQuantityFactory.
 */

class RISKMGR_DLL RiskAxis: public virtual IObject {
public:

    static CClassConstSP const TYPE;

    RiskAxis();
    ~RiskAxis();

    /**
     * "Un-frozen" IRiskAxis.
     */

    virtual IRiskAxisConstSP thawed() const = 0;

    /**
     * Register a scalar or vector IRiskAxis class for use with frozenSimple().
     *
     * The class must have exactly two (non-transient) fields: one of type
     * OutputName, the other either of type Void (for scalar risk axes) or of
     * type ExpiryWindow (for "vector" term-structured ones).
     *
     * For scalar risk axes which don't define a Void qualifier field, use
     * registerSimple(string, string), below.
     *
     * Example: PropertyRiskAxis<PROPERTY>::load().
     *
     * @param className            The name of the IRiskAxis class.  (Why not
     *                             the CClass itself?  Because this gets called
     *                             from the class's "load" function.)
     *                             
     * @param marketNameFieldName  The name of the field on the class which
     *                             defines the market data name with respect
     *                             to which the IRiskAxis is defined.
     *                             
     * @param qualifierFieldName   The name of the field on the class which
     *                             defines the qualifier with respect
     *                             to which the IRiskAxis is defined.
     */

    static void registerSimple(string className,
                               string marketNameFieldName,
                               string qualifierFieldName);

    /**
     * Register a "scalar" IRiskAxis class for use with frozenSimple().
     *
     * The class must have exactly one (non-transient) field, of type
     * OutputName.  For "vector" (term-structured) risk axes, use
     * registerSimple(string, string, string), above.
     *
     * @param className            The name of the IRiskAxis class.  (Why not
     *                             the CClass itself?  Because this gets called
     *                             from the class's "load" function.)
     *                             
     * @param marketNameFieldName  The name of the field on the class which
     *                             defines the market data name with respect
     *                             to which the IRiskAxis is defined.
     */

    static void registerSimple(string className,
                               string marketNameFieldName);

    /**
     * A suitable frozen representation for a simple IRiskAxis which has been
     * registered using registerSimple().
     *
     * Provides a default implementation of IRiskAxis::frozen() /
     * RiskAxis::thawed() for straightforward risk axes which can be defined
     * just by a market data name and perhaps an expiry.
     *
     * Returns either ScalarRiskAxis or ExpiryRiskAxis, both defined in
     * RiskAxis.cpp.
     *
     * Example: PropertyRiskAxis<PROPERTY>::frozen().
     */

    static RiskAxisConstSP frozenSimple(const IRiskAxis* it);
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
