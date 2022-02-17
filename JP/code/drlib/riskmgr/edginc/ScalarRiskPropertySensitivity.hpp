/**
 * @file ScalarRiskPropertySensitivity.hpp
 */

#ifndef QLIB_ScalarRiskPropertySensitivity_H
#define QLIB_ScalarRiskPropertySensitivity_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/Perturbation.hpp"
#include "edginc/IRiskProperty.hpp"
#include "edginc/ScalarPerNameShift.hpp"
#include "edginc/PerNameRiskPropertySensitivity.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(IResultsFunction)
FORWARD_DECLARE(ScalarRiskPropertySensitivity)

/**
 * A greek calculated by tweaking a non-term-structured "risk property"
 * (like spot) on each of the "names" in the market in turn.
 *
 * Semantics are described in the constructor ScalarRiskPropertySensitivity().
 *
 * See IRiskQuantityFactory for an overview of the "declarative" sensitivities
 * framework of which this class is a part, and RiskPropertySensitivity for
 * information about closely related classes.  (It's the equivalent of the
 * "old" ScalarShift.)
 */

class RISKMGR_DLL ScalarRiskPropertySensitivity:
       public PerNameRiskPropertySensitivity<Void>,
       public virtual IPerturbation,
       public virtual IScalarPerNameShift {

    static void load(CClassSP& clazz);

public:

    static CClassConstSP const TYPE;

protected:

    /**
     * Constructor
     *
     * This is really the first half of the constructor, the second half being
     * RiskPropertySensitivity::deriv() --- it's called after the object's
     * reflection FIELDs have been filled in so that it can refer to their
     * values.  The fields of the Deriv object you return from deriv() are
     * listed here as if they were constructor arguments, so that you can see
     * everything in one place.
     *
     * @param type             Reflection type of final subclass
     *
     * @param outputName       Packet name under which quantities computed for
     *                         this greek will be stored in the Results
     *                         (e.g. "PAR_SPREAD_RHO_PARALLEL" for
     *                         ParSpreadRhoParallel.cpp)
     *
     * @param shiftSize        How big to make the @a property tweaks used to
     *                         estimate @a derivative
     *
     * @param derivand         The quantity whose derivative we want to
     *                         estimate: generally fair value (i.e.
     *                         IResultsFunction::price()) but it could be
     *                         NAKED_BOND_PRICE or anything
     * 
     * @param property         The ExpiryWindow-qualified IRiskProperty with
     *                         respect to which the @a derivative is to be estimated
     *                         (e.g. RiskProperty<ParSpreadParallel> for
     *                         ParSpreadRhoParallel.cpp)
     *
     * @param derivative       The type of derivative to estimate, e.g.
     *                         IScalarDerivative::oneSided(),
     *                         IScalarDerivative::twoSided(),
     *                         IScalarDerivative::second()
     *
     * @param sensitivityUnit  Scaling factor applied to computed derivatives
     *                         before they're placed in the Results.  Normally
     *                         it should just be 1.0, but people want to see
     *                         some sensitivities reported in terms of basis
     *                         points (0.0001).
     *
     * @return  When the sensitivity is evaluated (typically via
     *    RiskQuantityEvaluator::storeResults()), the following derivatives will
     *    be returned from nameRiskQuantities(), estimated, and placed in
     *    the Results:
     *
     *   -  for each "name" in the world which implements
     *      property.subjectInterface()
     *      (i.e. ITweakableWithRespectTo<PROPERTY>),
     *
     *      -  the specified derivative of derivand with respect to
     *         property at that expiry, estimated by a perturbation
     *         of magnitude @a shiftSize, scaled by @a sensitivityUnit
     *
     *      -  ... placed in packet @a outputName / market data name
     */

    ScalarRiskPropertySensitivity(CClassConstSP type,
                                  const string& outputName,
                                  double shiftSize);

    /**
     * Constructor, with packet name for second derivative
     *
     * You can specify two derivatives to be computed for each market name
     * (typically IScalarDerivative::twoSided() and
     * IScalarDerivative::second()).  You'll also need to supply @a derivative2
     * and @a sensitivityUnit2 in the RiskPropertySensitivity::Deriv you
     * return from RiskPropertySensitivity::deriv().
     */

    ScalarRiskPropertySensitivity(CClassConstSP type,
                                  const string& outputName,
                                  const string& outputName2,
                                  double shiftSize);

public:

    ~ScalarRiskPropertySensitivity();

    /**
     * IPerturbation implementation (compatibility)
     */

    bool findAndShift(IObjectSP objectToShift, OutputNameConstSP name);

    /**
     * IScalarPerNameShift implementation (compatibility, for ImpliedScalarShift)
     */

    //@{

    IHypothesis::AlternateWorldSP appliedTo(OutputNameConstSP name,
                                            double shiftSize,
                                            IObjectSP world);
    double getShiftSize() const;

    OutputNameArrayConstSP allNames(const IObject* object) const;

    //@}
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
