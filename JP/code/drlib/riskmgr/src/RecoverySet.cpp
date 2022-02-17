//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RecoverySet.cpp
//
//   Description : Scenario for adjusting CDS Recovery
//
//   Author      : Mark A Robson
//
//   Date        : 26 July 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/RecoveryShiftBase.hpp"

DRLIB_BEGIN_NAMESPACE

class RecoverySet: public RecoveryShiftBase{
public:
    static CClassConstSP const TYPE;
    const static string NAME;

    /** check the inputs */
    virtual void validatePop2Object(){
        double newRecovery = getShiftSize();
        if (newRecovery < 0.0 || newRecovery > 1.0){
            throw ModelException("RecoverySet::validatePop2Object",
                                 "recovery must lie in [0,1]");
        }
    }

    /** Returns the new recovery level given the original one */
    virtual double applyShift(double unadjRecovery){
        setInitialValue(unadjRecovery);
        return getShiftSize();
    }

    /** Returns the original recovery level given the adjusted one */
    virtual double undoShift(double adjRecovery){
        return getInitialValue();
    }

    /** overridden: pretty meaningless for a scenario */
    virtual double divisor() const{
        return 1.0;
    }

    /** overridden: fails */
    virtual void calculate(TweakGroup* tweakGroup,
                           CResults*    results) {
        throw ModelException("RecoverySet::calculate", "Not supported");
    }
    
private:

    /** for reflection */
    RecoverySet() : RecoveryShiftBase(TYPE, NAME) 
    {}

    static IObject* defaultConstructor(){
        return new RecoverySet();
    }

    /** Invoked when SmallClass is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(RecoverySet, clazz);
        SUPERCLASS(RecoveryShiftBase);
        EMPTY_SHELL_METHOD(defaultConstructor);
    }
};


CClassConstSP const RecoverySet::TYPE = CClass::registerClassLoadMethod(
    "RecoverySet", typeid(RecoverySet), load);

const string RecoverySet::NAME = "RECOVERY_SET";

/** Included in RiskMgrLib::linkInClasses() to force link to include this */
bool RecoverySetLinkIn() {
    return RecoverySet::TYPE != NULL;
}


DRLIB_END_NAMESPACE
