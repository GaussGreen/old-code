//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Author      : Mark A Robson
//
//   Date        : 23 Dec 2005
//
//----------------------------------------------------------------------------

#ifndef QR_CCMLOSSCALCULATORBASE_HPP
#define QR_CCMLOSSCALCULATORBASE_HPP

#include "edginc/CreditMetricsLossCalculatorBase.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(CcmOnlyParameters);
FORWARD_DECLARE(SingleCreditAsset);

/** Base class to ease implementations of ITrancheLossCalculator for CCM */
class CONVOLUTION_DLL CCMLossCalculatorBase:
    public CreditMetricsLossCalculatorBase
{
public:
    virtual ~CCMLossCalculatorBase();

    /** Constructor - takes in full timeline to allow for optimisations. If
        computeCondCurve if false then lossConditional in loss() below will be
        set to zero. */
    CCMLossCalculatorBase(
        const DateTimeArray&           timeline, /* (I) */
        CreditTrancheLossConfigConstSP tranche,  /* (I) */
        CounterPartyCreditConstSP      cpty);    /* (I) */

    /** Retrieves CCM stlye engine params from Credit asset. Fails if they are
        not there or are of wrong type */
    static CcmOnlyParametersConstSP getCCMEngineParameters(
        SingleCreditAssetConstSP asset);


protected:
    DoubleArrayArray floorSurvivalProb;
};

DRLIB_END_NAMESPACE
#endif
