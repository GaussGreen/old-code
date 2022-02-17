/**
 * @file BasketDelta.cpp
 */

#ifndef DRLIB_BasketDelta_H
#define DRLIB_BasketDelta_H

#include "edginc/config.hpp"
#include "edginc/BasketSpot.hpp"
#include "edginc/ScalarShiftTwoSided.hpp"
#include "edginc/AssetSpotGreek.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * As Delta but for baskets.
 *
 * BasketDelta actually gets invoked whenever Delta is invoked, but right now
 * I'm documenting the new stuff not the old stuff ;).
 */

class RISKMGR_DLL BasketDelta: public ScalarShiftTwoSided,
                   public virtual Additive,
                   public virtual IAssetSpotGreek {

    static void load(CClassSP& clazz);

public:

    static CClassConstSP const TYPE;
    static const double DEFAULT_SHIFT;

private:

    bool useDeltaNames;

public:

    BasketDelta(double shiftSize);

    BasketDelta(double shiftSize,
                IModel* model, CControl* control, bool useDeltaNames);

    ~BasketDelta();

    const string& getSecondOrderSensOutputName() const;

    double divisor() const;
};

DECLARE(BasketDelta)

DRLIB_END_NAMESPACE

#endif
