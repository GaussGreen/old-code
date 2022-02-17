/**
 * @file RiskQuantityEvaluator.hpp
 */

#ifndef QLIB_RiskQuantityEvaluator_H
#define QLIB_RiskQuantityEvaluator_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/Results_forward.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(MultiTweakGroup)
FORWARD_DECLARE(RiskQuantityEvaluator)
FORWARD_DECLARE(RiskMapping)
FORWARD_DECLARE(NamedRiskQuantity)
FORWARD_DECLARE(IRiskQuantityFactory)

/**
 * Evaluates greeks
 *
 * This the class which does the work of tweaking/pricing/calculating to
 * evaluate greeks written in the "declarative" sensitivities framework (see
 * IRiskQuantityFactory for an overview).
 *
 * There are two entry points:
 *
 *   -   values() interrogates a list of greeks for the RiskQuantity's they
 *       want evaluated, evaluates them, and returns them along with their
 *       values
 *
 *   -   storeResults() goes a step further and populates a Results array
 *       with the RiskQuantity values
 */

class RISKMGR_DLL RiskQuantityEvaluator: public CObject {

    static void load(CClassSP& clazz);
    RiskQuantityEvaluator(const RiskQuantityEvaluator& rhs);
    RiskQuantityEvaluator& operator=(const RiskQuantityEvaluator& rhs);

public:

    static CClassConstSP const TYPE;

private:

    bool riskMappingEnabled;
    RiskMappingConstSP riskMapping;

public:

    /**
     * Constructor
     *
     * If @a riskMappingEnabled, and getMarket() finds suitable
     * RiskMappingMatrix's in the MarketData, then they will be used to
     * evaluate RiskQuantity's specified in terms of Black-Scholes
     * tweaks using "equivalent" IDynamicsParameter tweaks instead.
     */

    RiskQuantityEvaluator(bool riskMappingEnabled = false);
	RiskQuantityEvaluator(RiskMappingSP riskMapping);
    ~RiskQuantityEvaluator();

    /**
     * Retrieve necessary data from the market
     *
     * Called at initialisation time from Control::getMarket().
     * It constructs a RiskMapping, if riskMappingEnabled.
     */

    void getMarket(IModelConstSP model, MarketDataConstSP market,
                   IInstrumentCollectionConstSP instruments);

    FORWARD_DECLARE(Value)

    /**
     * Value computed for a NamedRiskQuantity for a certain instrument, or
     * information about what went wrong calculating it
     *
     * You only need to know about this class if you call
     * RiskQuantityEvaluator::values() ---
     * RiskQuantityEvaluator::storeResults() puts the values straight into
     * Results dictionaries.
     *
     * The important things here are instrument, riskQuantity and value(): the
     * former tell you what was evaluated, the latter its value.  If something
     * went wrong oops tells you what.  For diagnostic purposes you also have
     * access to the intermediate vals and dsts used in the calculation.
     */

    class RISKMGR_DLL Value: public CObject {

        Value(const Value&);
        Value &operator =(const Value&);

        static void load(CClassSP& clazz);

        mutable double _value; // $unregistered
        mutable bool _valueValid; // $unregistered

    public:

        static CClassConstSP const TYPE;

        /**
         * The RiskQuantity for which we give a value
         */

        NamedRiskQuantitySP riskQuantity;

        /**
         * Number of instrument within IInstrumentCollection for which we're
         * evaluating riskQuantity
         */

        int instrument;

        /**
         * The values of the HypotheticalQuantity's needed for computing
         * riskQuantity.
         *
         * Specifically,
         *
         *    -  vals[i] gets set to the value of the i'th of riskQuantity's
         *       RiskQuantity::parameters(), by applying its
         *       HypotheticalQuantity::hypothesis() to the world and then
         *       calling its HypotheticalQuantity::quantity() function
         *
         *    -  dsts[i] gets set to the distance of that i'th world from
         *       the base case world (i.e. appropriate "divisor" but without
         *       scaling factor)
         */

        DoubleArraySP vals, dsts;

        /**
         * Any problems encountered in putting the world into an appropriate
         * state, pricing, or computing hypothetical quantities
         */

        UntweakableConstSP oops;

        /**
         * Whether the hypotheses involved in calculating this value all made
         * sense
         *
         * If false, some of the tweaks applied in calculating the value were
         * "not applicable" because e.g. no objects were found in the world
         * which fell into their domain.
         */

        bool hypothesesApplicable;

        /**
         * Whether the instrument involved in calculating this value all made
         * sense
         *
         * If false, some quantity involved in calculating the value is not
         * defined for the instrument concerned (for instance an output
         * request which it doesn't support)
         */

        bool instrumentApplicable;

        /**
         * Whether the calculation of the quantity involved RiskMapping
         */

        bool usedRiskMapping;

        /**
         * Constructor.
         */

        Value(NamedRiskQuantitySP riskQuantity, int instrument);

        /**
         * Constructing returning a smartPtr.
         */

        static ValueSP SP(NamedRiskQuantitySP riskQuantity, int instrument);

        ~Value();

        /**
         * Value obtained for riskQuantity
         *
         * Throws an exception if something went wrong during the calculation
         * (untweakable/not applicable/zero divisor ...).
         */

        double value();

        /**
         * False if value() is going to throw an exception
         */

        bool ok();

        /**
         * Evaluate the riskQuantity and put its value in the appropriate
         * Results
         *
         * Roughly:
         *
         * -     resultss[instrument][riskQuantity->resultsName] <- value())
         */

        void storeResult(CResultsArraySP resultss);
    };

    /**
     * The RiskQuantity's generated by @a greeks, together with their values
     * computed against @a world
     *
     * This is the top-level entry point for calculating sensitivities in the
     * declarative framework.
     *
     * We take a MultiTweakGroup (= Model and an IInstrumentCollection to
     * price with it), and a list of greeks, say Delta and VegaPointwise.
     *
     * We interrogate the greeks for the NamedRiskQuantity's they want
     * to evaluate, using see IRiskQuantityFactory::riskQuantity().
     *
     * Then we collate the IHypothesis's mentioned in the RiskQuantity's, and
     * run through them generating alternate worlds, recording their distances
     * from the base world, and evaluating instrument prices (or whatever the
     * RiskQuantity's in question want).
     *
     * We write the hypothetical prices and corresponding distances for each
     * RiskQuantity into a Value record, and return them all: for instance,
     *
     *    -  DELTA IBM for instrument 0 = 16
     *    -  DELTA IBM for instrument 1 = -3
     *    -  VEGA_POINTWISE IBM 3M for instrument 0 = 22
     *    -  VEGA_POINTWISE IBM 6M for instrument 0 = 0.5
     *    -  VEGA_POINTWISE INTC 3M for instrument 0 = 1.2
     *    -  ...
     */

    ValueArraySP values(
        IRiskQuantityFactoryArrayConstSP greeks,
        MultiTweakGroupSP world) const;

    /**
     * Evaluate the RiskQuantity's generated by @a greeks against @a world, and
     * put them in @a resultss
     *
     * As values(), but writes the computed RiskQuantity's directly into resultss
     * rather than returning them as Value's.  resultss[i] gets populated with
     * the numbers relating to the i'th instrument in the IInstrumentCollection
     * world->getInstruments().
     *
     * Called from Control::calculateMulti().
     */

    void storeResults(IRiskQuantityFactoryArrayConstSP greeks,
                      MultiTweakGroupSP world,
                      CControlSP control,
                      CResultsArraySP resultss) const;

    /**
     * Evaluate the RiskQuantity's generated by @a greeks against @a world but
     * only for the purpose of computing exposure characteristics, and
     * put them in @a resultss b
     *
     * As values(), but writes the computed RiskQuantity's directly into resultss
     * rather than returning them as Value's.  resultss[i] gets populated with
     * the numbers relating to the i'th instrument in the IInstrumentCollection
     * world->getInstruments().
     *
     * Called from Control::calculateMultiExposures().
     */

    void storeExposureResults(IRiskQuantityFactoryArrayConstSP greeks,
        MultiTweakGroupSP world,
        CControlSP control,
        CResultsArraySP resultss) const;
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
