//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RecoveryTweak.hpp
//
//   Description : Sensitivity to Recovery
//
//   Date        : 19 October 2005
//
//----------------------------------------------------------------------------

#ifndef QLIB_RECOVERYTWEAK_HPP
#define QLIB_RECOVERYTWEAK_HPP

#include "edginc/RecoveryShiftBase.hpp"
#include "edginc/Class.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE

//// Sensitivity to CDS Recovery
class RISKMGR_DLL RecoveryTweak: public RecoveryShiftBase{
public:
    static CClassConstSP const TYPE;
    const static string NAME;
    const static double DEFAULT_SHIFT;

    /** constructor with explicit shift size */
    RecoveryTweak(double shiftSize);

    /** Returns the new recovery level given the original one */
    virtual double applyShift(double unadjRecovery);

    /** Returns the original recovery level given the adjusted one */
    virtual double undoShift(double adjRecovery);

private:
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class RISKMGR_DLL Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar 
    {
    public:
        virtual Sensitivity* createDefault();
        virtual Sensitivity* createScalar(double shiftSize);
    };

    /** for reflection */
    RecoveryTweak();

    static IObject* defaultRecoveryTweak();

    /** Invoked when SmallClass is 'loaded' */
    static void load(CClassSP& clazz);
};

DRLIB_END_NAMESPACE

#endif
