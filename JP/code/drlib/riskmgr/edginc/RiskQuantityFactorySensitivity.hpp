/**
 * @file RiskQuantityFactorySensitivity.hpp
 */

#ifndef QLIB_RiskQuantityFactorySensitivity_H
#define QLIB_RiskQuantityFactorySensitivity_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/Sensitivity.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/IRiskQuantityFactory.hpp"
#include "edginc/ICompatibilitySensitivity.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(RiskQuantityFactorySensitivity)

/**
 * A Sensitivity implemented using the "declarative" risk management scheme
 *
 * See IRiskQuantityFactory for an overview of the framework; this class
 * provides an implementation of Sensitivity in terms of IRiskQuantityFactory,
 * so that sensitivities can be redefined in the "new" way and still be
 * used unchanged in existing contexts.
 *
 * The key method implemented here is calculate() (= Sensitivity::calculate()).
 * Note that it hardly ever gets called under normal circumstances: the method
 * usually responsible for triggering calculation of sensitivities,
 * Control::calculateMulti(), treats RiskQuantityFactorySensitivity's
 * separately from other Sensitivity's, evaluating them all together using
 * RiskQuantityEvaluator::storeResults() rather than calculate() in order to exploit
 * its greater efficiency in the multi-instrument case.
 *
 * The only abstract method which subclasses need to provide is
 * nameRiskQuantities().
 *
 * Interesting subclasses are RiskPropertySensitivity, which is a sensitivity
 * to a "risk property" (like spot or vol) of names in the market,
 * PerNameRiskPropertySensitivity, and concretely VegaPointwise,
 * ParSpreadRhoParallel.cpp etc.
 */

class RISKMGR_DLL RiskQuantityFactorySensitivity:
        public Sensitivity,
        public virtual IRiskQuantityFactory,
        public virtual ICompatibilitySensitivity {

    static void load(CClassSP& clazz);

public:

    static CClassConstSP const TYPE;

protected:

    const StringArrayConstSP outputNames;

    /**
     * ICompatibilitySensitivity implementation
     *
     * These methods support an interim mechanism which enables
     * instrument/model logic based on interrogating
     * Control::getCurrentSensitivity() to work more or less unchanged.
     * Those mechanisms will likely be overhauled to work directly
     * (and more robustly) in terms of the difference between the
     * world being priced and the base-base world.
     */

    //@{

    void setCurrentHypothesis(AtomicHypothesisArraySP history,
                              AbstractPropertyTweakHypothesisConstSP current,
                              double oldValue);

    mutable OutputNameConstSP marketDataName; // $unregistered

    //@}

    /**
     * Constructor
     *
     * @param type             Reflection type of final subclass
     *
     * @param outputName       Packet name under which quantities computed for
     *                         this greek will be stored in the Results
     *                         (e.g. "PAR_SPREAD_RHO_PARALLEL" for
     *                         ParSpreadRhoParallel.cpp)
     *
     * @param outputName2      Packet name under which quantities computed for
     *                         the optional second derivative will be stored
     *                         in the results --- see the alternative
     *                         constructor on
     *                         RiskPropertySensitivity<QUALIFIER>::Deriv
     */

    RiskQuantityFactorySensitivity(CClassConstSP type,
                                   const string& outputName,
                                   const string& outputName2 = "");

    /**
     * The RiskQuantity's which should be computed for a given base-state
     * world
     *
     * Each NamedRiskQuantity you return designates an individual scalar to be
     * computed as an output of this sensitivity, like "IBM 6M vega", "IBM 1Y
     * vega", "JPM 1M vega".  Note that you don't have to worry about things
     * like returning "not applicable" in case there are no other outputs,
     * because that's handled by the outer riskQuantities() entry point.
     *
     * The canonical implementation is
     * PerNameRiskPropertySensitivity::nameRiskQuantities().
     */

    virtual NamedRiskQuantityArraySP nameRiskQuantities(
        MultiTweakGroupConstSP world, RiskMappingConstSP riskMapping) const = 0;

    /**
     * Implementation of Sensitivity::calculate() in terms of IRiskQuantityFactory
     *
     * As noted in the class doc above, this method allows "declarative"
     * sensitivities behave as drop-in replacements for Sensitivity; but it's
     * hardly used, because Control::calculateMulti() bypasses it and goes
     * straight to RiskQuantityEvaluator::storeResults(), in order to exploit its
     * greater efficiency in the multi-instrument case.
     */

    virtual void calculate(TweakGroup* tweakGroup, Results* results);
    friend class FlexibleSensitivity; // aaaagh

public:

    ~RiskQuantityFactorySensitivity();

    /**
     * The packet name against which the numbers generated by this sensitivity
     * should be stored in the Results ("DELTA")
     *
     * Comes straight from the @a outputName argument to the
     * RiskQuantityFactorySensitivity() constructor.
     */

    virtual const string& getSensOutputName() const;

    /**
     * Implementation of Sensitivity::scaleResult()
     */

    virtual void scaleResult(Results* results, double scaleFactor) const;

    /**
     * Implementation of Sensitivity::addResult()
     */

    virtual void addResult(Results* results, const Results* resultsToAdd,
                           double scaleFactor) const;

    /**
     * The RiskQuantity's which should be computed for a given base-state
     * world
     *
     * Returns a list of all the individual scalar quantities (like "IBM 6M
     * vega", "IBM delta", "JPM 1M vega") which should be computed.
     *
     * These come from the abstract nameRiskQuantities() method but if the list is
     * empty we insert a single "not applicable".
     */

    virtual NamedRiskQuantityArraySP riskQuantities(
        MultiTweakGroupConstSP world, RiskMappingConstSP riskMapping) const;

    /**
     * Mechanism for handling the case in which the set of desired
     * RiskQuantity's can't be determined from the base-state world
     *
     * In all current cases it can, so the default implementation just returns
     * a null smartPtr.  See IRiskQuantityFactory::lazies().
     */

    virtual LazyRiskQuantityFactoryArraySP lazies(MultiTweakGroupConstSP world) const;

    /**
     * An RiskQuantityFactorySensitivity which generates a single NamedRiskQuantity
     *
     * This is used for generating error messages in the results:
     * see for instance FieldSensitivityDefinition::sensitivity().
     */

    static RiskQuantityFactorySensitivitySP singleton(
        const string& packetName,
        NamedRiskQuantityConstSP rq);
};

DRLIB_END_NAMESPACE

#endif
