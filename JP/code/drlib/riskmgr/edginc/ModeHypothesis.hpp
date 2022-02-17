/**
 * @file ModeHypothesis.hpp
 */

#ifndef QLIB_ModeHypothesis_H
#define QLIB_ModeHypothesis_H

#include "edginc/AtomicHypothesis.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(ModeHypothesis)
FORWARD_DECLARE(MultiTweakGroup)

/**
 * A function which puts an instrument/model assembly into a different pricing
 * "mode".
 *
 * Some Greeks (e.g. LiquiditySpreadRhoPointwise) are defined with respect to
 * prices computed once the model has been put into an alternate "mode".  This
 * interface allows easy definition of such a hypothesis.  See
 * LiquiditySpreadRhoPointwise.cpp for an example.
 *
 * You need to provide applyTo(MultiTweakGroupSP) and undo(MultiTweakGroupSP).
 * Recall that MultiTweakGroup is the multi-instrument version of
 * TweakGroup.
 */

class RISKMGR_DLL ModeHypothesis: public AtomicHypothesis {

    ModeHypothesis();
    static IObject* defaultOne();
    static void load(CClassSP& clazz);

    friend class ModeHypothesis_AltWorld;

public:

    static CClassConstSP const TYPE;

protected:

    virtual double applyToWorld(MultiTweakGroupSP world,
                                bool* changed) const = 0;
    virtual void undo(MultiTweakGroupSP world) const = 0;

public:

    IRiskAxisConstSP axis() const;

    double axisCoefficient() const;

    double applyTo(IObjectSP world, bool* changed = 0) const;

    IHypothesis::AlternateWorldSP appliedTo(IObjectSP world) const;

    ModeHypothesis(CClassConstSP type);
    ~ModeHypothesis();
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
