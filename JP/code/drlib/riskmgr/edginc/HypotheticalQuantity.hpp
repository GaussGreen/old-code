/**
 * @file HypotheticalQuantity.hpp
 */

#ifndef DRLIB_HypotheticalQuantity_H
#define DRLIB_HypotheticalQuantity_H

#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(IHypothesis)
FORWARD_DECLARE(IResultsFunction)
FORWARD_DECLARE(HypotheticalQuantity)

/**
 * A quantity (typically instrument fair value) computed in a hypothetical
 * alternate state of the world (= instrument/model assembly).
 *
 * See RiskQuantity for a description of what this class is for, and
 * IRiskQuantityFactory for an overview of the "declarative" sensitivities
 * framework of which it forms a part.
 *
 * Note that while HypotheticalQuantity describes what hypothesis() to apply to
 * the world, and the quantity() to evaluate in the resulting state, it doesn't
 * actually have a method to do that: the work is performed by
 * RiskQuantityEvaluator::storeResults(), so that it can be shared across the many
 * HypotheticalQuantity's arising from an array of IRiskQuantityFactory's.
 */

class RISKMGR_DLL HypotheticalQuantity: public CObject {

    HypotheticalQuantity(const HypotheticalQuantity& rhs);
    HypotheticalQuantity& operator=(const HypotheticalQuantity& rhs);
    static void load(CClassSP& clazz);

public:

    static CClassConstSP const TYPE;

private:

    IHypothesisConstSP _hypothesis;
    IResultsFunctionConstSP _quantity;

public:

    /**
     * Constructor.
     *
     * @param hypothesis  IHypothesis defining the alternate state into which
     *                    the world is to be put before evaluating @a quantity
     *
     * @param quantity    IResultsFunction to be applied to the pricing Results
     *                    obtained in the alternate state of the world defined
     *                    by @a hypothesis
     */

    HypotheticalQuantity(IHypothesisConstSP hypothesis,
                         IResultsFunctionConstSP quantity);

    /**
     * Constructor returning a smartPtr.
     */

    static HypotheticalQuantitySP SP(IHypothesisConstSP hypothesis,
                                     IResultsFunctionConstSP quantity);

    ~HypotheticalQuantity();

    /**
     * IHypothesis defining the alternate state into which the world is to be
     * put before evaluating quantity()
     */

    IHypothesisConstSP hypothesis() const;

    /**
     * IResultsFunction to be applied to the pricing Results obtained in the
     * alternate state of the world defined by hypothesis()
     *
     * It's expressed in the form of an IResultsFunction function for
     * extracting a number from a Results dictionary after pricing the
     * instrument.  Typically it's IResultsFunction::price(), which retrieves
     * "fair value", but it could for instance be a NAKED_BOND_PRICE output
     * request, or whatever.
     */

    IResultsFunctionConstSP quantity() const;
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
