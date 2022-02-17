/**
 * @file Delta.hpp
 */

#ifndef DRLIB_Delta_H
#define DRLIB_Delta_H

#include "edginc/Spot.hpp"
#include "edginc/ScalarShiftTwoSided.hpp"
#include "edginc/AssetSpotGreek.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(Delta)

/**
 * A greek calculated by tweaking spot on all or some of the "names" in the
 * market which have it.
 *
 * Delta has been largely rewritten to use the "declarative" sensitivities
 * framework: see IRiskQuantityFactory for an overview.  For the moment,
 * however, it and its sister BasketDelta are <I>not</I> implemented using
 * ScalarRiskPropertySensitivity as they really ought to be.  This is because
 * there are some fiddly mechanisms (like MonteCarlo "quick greeks") which are
 * written in terms of ScalarShift and I don't want to have to fix them up
 * before committing this stuff.  Also the mechanism whereby Delta
 * automatically invokes BasketDelta is a little specialised---maybe we
 * will rethink it at some point.
 *
 * The key method here is riskQuantities()
 * (IRiskQuantityFactory::riskQuantities()).  calculate()
 * (Sensitivity::calculate()) is implemented, but the central point where
 * sensitivity calculations are triggered, Control::calculateMulti(), doesn't
 * use it.
 */

class RISKMGR_DLL Delta: public ScalarShiftTwoSided,
             public virtual Additive,
             public virtual IAssetSpotGreek {

    static void load(CClassSP& clazz);

public:

    static CClassConstSP const TYPE;

    static const string NAME;
    static const string SECOND_ORDER_NAME;
    static const double DEFAULT_SHIFT;
    static const double MINIMUM_SHIFT;

    /**
     * Constructor.
     *
     * @param shiftSize  Relative perturbation to be made to spot
     *                   when estimating derivative
     *                   (spot' = spot * (1 +shiftSize))
     */

    Delta(double shiftSize);
    ~Delta();

    /**
     * Constructor (deprecated).
     *
     * Used to implement AlterControl.
     */

    Delta(double shiftSize, IModel* model, CControl* control);

    /**
     * Implementation of Sensitivity::calculate()
     */

    void calculate(TweakGroup* tweakGroup, CResults* results);

    /**
     * Packet name under which to store second deriv values ("GAMMA")
     */

    const string& getSecondOrderSensOutputName() const;

    /**
     * Divisor to use in calculating deriv =~ (f(x') - f(x)) / divisor
     *
     * This is the absolute amount by which spot was tweaked (not the relative
     * perturbation specified in the Delta() constructor).  Hence the
     * derivative is reported in absolute not relative units.
     */

    double divisor() const;

    NamedRiskQuantityArraySP riskQuantities(
        MultiTweakGroupConstSP world, RiskMappingConstSP riskMapping) const;
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
