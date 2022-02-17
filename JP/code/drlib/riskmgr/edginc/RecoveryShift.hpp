//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RecoveryShift.hpp
//
//   Description : Scenario for adjusting CDS Recovery
//
//----------------------------------------------------------------------------

#ifndef QLIB_RECOVERYSHIFT_HPP
#define QLIB_RECOVERYSHIFT_HPP

#include "edginc/Class.hpp"
#include "edginc/Maths.hpp"
#include "edginc/RecoveryShiftBase.hpp"

DRLIB_BEGIN_NAMESPACE

class RISKMGR_DLL RecoveryShift: public RecoveryShiftBase {
public:
    static CClassConstSP const TYPE;
    const static string NAME;

    virtual ~RecoveryShift();
    // Explicit constructor
    RecoveryShift(double shiftSize, bool relativeShift);

    /** check the inputs */
    virtual void validatePop2Object();

    /** Returns the new recovery level given the original one */
    virtual double applyShift(double unadjRecovery);

    /** Returns the original recovery level given the adjusted one */
    virtual double undoShift(double adjRecovery);

    /** overridden: pretty meaningless for a scenario */
    virtual double divisor() const;

    /** overridden: fails */
    virtual void calculate(TweakGroup* tweakGroup,
                           CResults*    results);
    
protected:
    RecoveryShift(CClassConstSP clazz, const string& outputName);

private:
    bool   relativeShift; // optional, default false
    double recoveryFloor; // optional, default 0.0
    double recoveryCap;   // optional, default 1.0

    /** for reflection */
    RecoveryShift();
    static IObject* defaultConstructor();
    /** Invoked when SmallClass is 'loaded' */
    static void load(CClassSP& clazz);
};

DRLIB_END_NAMESPACE

#endif
