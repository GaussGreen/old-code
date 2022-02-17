//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : RecoveryTweakWithShift.hpp
//
//   Description : Applies two tweaks to shift/set the Recovery Rate:
//                 First shifts the recovery rate and then calculates the
//                 RecoveryTweak greek
//
//   Author      : Jose Hilera
//
//   Date        : 19 October 2005
//
//----------------------------------------------------------------------------

#ifndef QLIB_RECOVERYDOUBLETWEAK_HPP
#define QLIB_RECOVERYDOUBLETWEAK_HPP

#include "edginc/Class.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/RecoveryShiftBase.hpp"

DRLIB_BEGIN_NAMESPACE


class RISKMGR_DLL RecoveryTweakWithShift : public RecoveryShiftBase {
public:
    static CClassConstSP const TYPE;
    const static string NAME;

    virtual ~RecoveryTweakWithShift();
    // Explicit constructor
    RecoveryTweakWithShift(double shiftSize, bool relativeShift);

    /** Returns the new recovery level given the original one */
    virtual double applyShift(double unadjRecovery);

    /** Returns the original recovery level given the adjusted one */
    virtual double undoShift(double adjRecovery);

    virtual void calculate(TweakGroup* tweakGroup, CResults* results);

private:
    double initialShiftSize; // Size of the initial shift
    bool   relativeShift;  // optional, default false

    // for reflection 
    RecoveryTweakWithShift();
    static IObject* defaultConstructor();
    /** Invoked when SmallClass is 'loaded' */
    static void load(CClassSP& clazz);
};

DRLIB_END_NAMESPACE

#endif
