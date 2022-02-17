//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : LongTermSqueezeTweak.cpp
//
//   Description : LONG_TERM_SQUEEZE_TWEAK sensitivity
//
//   Author      : Oliver Brockhaus
//
//   Date        : 05 Mar 2003
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/LongTermSqueezeTweak.hpp"
#include "edginc/Maths.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE

LongTermSqueezeTweak::Shift::~Shift(){} // empty
LongTermSqueezeTweak::RestorableShift::~RestorableShift(){} // empty

/** Sens Control for LONG_TERM_SQUEEZE_TWEAK       */
const string LongTermSqueezeTweak::NAME = "LONG_TERM_SQUEEZE_TWEAK";
const double LongTermSqueezeTweak::DEFAULT_SHIFT = 0.01;

/** constructor with explicit shift size */
LongTermSqueezeTweak::LongTermSqueezeTweak(double shiftSize):
    ScalarShift(TYPE, NAME, shiftSize), useShiftSizeSign(false) {
    validatePop2Object();
}
LongTermSqueezeTweak::LongTermSqueezeTweak(double shiftSize, bool useShiftSizeSign):
    ScalarShift(TYPE, NAME, shiftSize), useShiftSizeSign(useShiftSizeSign){
    validatePop2Object();
}

/** for reflection */
LongTermSqueezeTweak::LongTermSqueezeTweak(): 
    ScalarShift(TYPE, NAME, DEFAULT_SHIFT), useShiftSizeSign(false){
}

/** Validation */
void LongTermSqueezeTweak::validatePop2Object(){
    static const string method("LongTermSqueezeTweak::validatePop2Object");
    double shiftSize = getShiftSize();
    if ( (shiftSize > 1.0) || (shiftSize < -1.0) ){
        throw ModelException( method, "LongTermSqueezeTweak "+Format::toString(shiftSize)+
                              "must be between -1.0 and 1.0.");
    }
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double LongTermSqueezeTweak::divisor() const{
    static const char routine[] = "LongTermSqueezeTweak::divisor";    
    double divisor;
    try{
        // find initial value of correlation squeeze - stored in sens control
        bool negate = false;
        if (!useShiftSizeSign) {
            try{
                double initVal = getInitialValue();
                if (!(Maths::isNegative(initVal))) {
                    negate = true;
                }
            } catch (ModelException& e){
                e.addMsg("Initial value of correlation squeeze has not been set"
                         " in the LONG_TERM_SQUEEZE_TWEAK calculation");
                throw e;
            }
        }
        // then just scale by shift size
        double shiftSize = getShiftSize();
        if (Maths::isZero(shiftSize)) {
            throw ModelException(routine, "Shift size is zero");
        }
        // report a 0.1 move in correlation sqeeze (in [-1,1]); shift towards 0.0
        divisor = 10.0 * shiftSize; 
        if (negate) {
            divisor = -divisor;
        }
    } catch (ModelException& e){
        e.addMsg(routine);
        throw e;
    }
    return divisor;
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP LongTermSqueezeTweak::shiftInterface() const{
    return Shift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP LongTermSqueezeTweak::restorableShiftInterface() const{
    return RestorableShift::TYPE;
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity.  The object must implement the
    VegaParallel.Shift interface */
bool LongTermSqueezeTweak::nameMatches(const OutputName&         name,
                        IObjectConstSP            obj){
    // cast obj to LongTermSqueezeTweak::Shift and then invoke name method
    const Shift& LongTermSqueezeTweakObj =
        dynamic_cast<const Shift&>(*obj);
    return name.equals(LongTermSqueezeTweakObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to
    this sensitivity to the supplied list */
void LongTermSqueezeTweak::appendName(OutputNameArray&          namesList,
                       IObjectConstSP            obj){
    // cast obj to LongTermSqueezeTweak::Shift and then invoke name method
    const Shift& LongTermSqueezeTweakObj =
        dynamic_cast<const Shift&>(*obj);
    OutputNameSP outputName(new OutputName(LongTermSqueezeTweakObj.sensName(this)));
    namesList.push_back(outputName);
}

bool LongTermSqueezeTweak::shift(IObjectSP obj) {
    // cast obj to LongTermSqueezeTweak::Shift and then invoke shift method
    Shift& LongTermSqueezeTweakObj =
        dynamic_cast<Shift&>(*obj);
    return LongTermSqueezeTweakObj.sensShift(this);
}

    
void LongTermSqueezeTweak::restore(IObjectSP obj) {
    // cast obj to LongTermSqueezeTweak::Shift and then invoke restore method
    RestorableShift& LongTermSqueezeTweakObj =
        dynamic_cast<RestorableShift&>(*obj);
    LongTermSqueezeTweakObj.sensRestore(this);
}

/** Returns true if implementations of Shift should override the sign of
    the shift size when making the tweak */
bool LongTermSqueezeTweak::overrideShiftSizeSign() const{
    return !useShiftSizeSign;
}

class LongTermSqueezeTweakHelper{
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new LongTermSqueezeTweak(LongTermSqueezeTweak::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new LongTermSqueezeTweak(shiftSize);
        }
    };
    
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(LongTermSqueezeTweak, clazz);
        SUPERCLASS(ScalarShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultLongTermSqueezeTweak);
        FIELD_NO_DESC(useShiftSizeSign);
        FIELD_MAKE_TRANSIENT(useShiftSizeSign);
        // register how to build our sensitivity
        SensitivityFactory::addSens(LongTermSqueezeTweak::NAME, 
                                    new Factory(), 
                                    new LongTermSqueezeTweak(LongTermSqueezeTweak::DEFAULT_SHIFT, false),
                                    LongTermSqueezeTweak::Shift::TYPE);
    }

    static IObject* defaultLongTermSqueezeTweak(){
        return new LongTermSqueezeTweak();
    }
};

CClassConstSP const LongTermSqueezeTweak::TYPE = CClass::registerClassLoadMethod(
    "LongTermSqueezeTweak", typeid(LongTermSqueezeTweak), LongTermSqueezeTweakHelper::load);

CClassConstSP const LongTermSqueezeTweak::Shift::TYPE =
CClass::registerInterfaceLoadMethod(
    "LongTermSqueezeTweak::Shift", typeid(LongTermSqueezeTweak::Shift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(LongTermSqueezeTweak::RestorableShift, clazz);
    EXTENDS(LongTermSqueezeTweak::Shift);
}
    
CClassConstSP const LongTermSqueezeTweak::RestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "LongTermSqueezeTweak::RestorableShift", typeid(LongTermSqueezeTweak::RestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE






