/**
 * @file Correl.hpp
 */

#ifndef DRLIB_Correl_H
#define DRLIB_Correl_H

#include "edginc/Void.hpp"
#include "edginc/IRestorableWithRespectTo.hpp"

DRLIB_BEGIN_NAMESPACE

struct Correl;
typedef smartConstPtr<Correl> CorrelConstSP;

/**
 * Tag class denoting "correlation" as a property of market names
 *
 * Currently, the only object which has this property (i.e. implements
 * ITweakableWithRespectTo<Correl> is Correlation; the sensitivities defined
 * with respect to this property are Phi, FXPhi, PhiPlus, CorrelationSqueeze,
 * and PhiParallel.
 *
 * See Spot for an account of what these tag classes are for.
 */

struct RISKMGR_DLL Correl: CObject {
    static CClassConstSP const TYPE;

    /**
     * How correlations can be tweaked
     *
     * Under ABSOLUTE_INWARD, shiftSize means the amount by which the correl is
     * moved towards zero.  Under ABSOLUTE, shiftSize is just the amount to
     * move the correl.  Under SQUEEZE, correl is moved by a proportion
     * shiftSize from where it is now towards 1.
     */

    enum Operation { ABSOLUTE_INWARD, ABSOLUTE, SQUEEZE };

    /**
     * Tag class denoting "correlation" as a property of market names
     *
     * These parameters get used in Correlation::sensShift() and
     * Correlation::sensName().
     *
     * @param asset   Whether the property encompasses asset-asset (eg EQ-EQ)
     *                correlations
     * @param fx      Whether the property encompasses FX-FX and FX-other
     *                correlations
     * @param other   Whether the property encompasses other (eg credit spreads-
     *                rates) correlations
     * @param op      How correlations are tweaked (see Operation)
     */

    Correl(bool asset = true, bool fx = true, bool other = true,
           Operation op = ABSOLUTE_INWARD);

    /**
     * Constructor returning smartPtr
     */

    static CorrelConstSP SP(bool asset, bool fx, bool other, Operation op);
    ~Correl();

    /**
     * Whether the property encompasses asset-asset (eg EQ-EQ) correlations
     *
     * Read by Correlation::sensShift().
     */

    bool asset;

    /**
     * Whether the property encompasses FX-FX and FX-other correlations
     *
     * Read by Correlation::sensShift().
     */

    bool fx;

    /**
     * Whether the property encompasses other (eg credit spreads-rates)
     * correlations
     *
     * Read by Correlation::sensShift().
     */

    bool other;

    /**
     * How the correlations are to be tweaked
     *
     * Read by Correlation::sensShift().
     */

    Operation op;

    /**
     * Info needed in addition to market data name to make a one-dimensional
     * "risk axis" for correlation (none)
     */

    typedef Void Qualifier;

    enum {
        /**
         * Correl is continuous, i.e. tweaks to it can be made arbitrarily small
         *
         * This is carried through to Sensitivity::discreteShift()
         * via RiskProperty<PROPERTY>::discrete().
         */

        discrete = 0
    };
};

#ifndef QLIB_CORREL_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL ITweakableWithRespectTo<Correl>);
EXTERN_TEMPLATE(class RISKMGR_DLL IRestorableWithRespectTo<Correl>);
#endif

DRLIB_END_NAMESPACE

#endif
