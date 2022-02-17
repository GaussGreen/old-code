/**
 * @file AtomicHypothesis.hpp
 */

#ifndef QLIB_AtomicHypothesis_H
#define QLIB_AtomicHypothesis_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/IHypothesis.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(AtomicHypothesis)

/**
 * An IHypothesis which isn't comprised of simpler component hypotheses.
 *
 * Almost all IHypothesis implementations are descended from this class, the
 * only exception being CompoundHypothesis.  This ensures that a
 * CompoundHypothesis is always a one-level list, never a multi-level tree
 * of compounds within compounds.
 *
 * See IRiskQuantityFactory for the big picture of the subsystem in which
 * this class is used.
 */

class RISKMGR_DLL AtomicHypothesis: public CObject,
                                    public virtual IHypothesis {

    static void load(CClassSP& clazz);

public:

    static CClassConstSP const TYPE;

    AtomicHypothesis(CClassConstSP type);
    ~AtomicHypothesis();

    /**
     * Number of atomic hypotheses to which this reduces: returns 1.  Don't
     * override it.
     */

    int numAtomics() const;

    /**
     * Atomic hypotheses to which this reduces: throws if i != 0 and returns
     * 'this'.  Don't override it.
     */

    AtomicHypothesisConstSP atomic(int i) const;

    IHypothesis::IDistanceMetricConstSP distanceMetric() const;

	/**
     * Second order approximation for risk mapping
     */

    int getApproxOrder() const;
	
	void setApproxOrder(int order);

private:

	int approxOrder;

};

DRLIB_END_NAMESPACE

#endif
