//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ShortTermSqueezeTweak.cpp
//
//   Description : SHORT_TERM_SQUEEZE_TWEAK sensitivity
//
//   Author      : Oliver Brockhaus
//
//   Date        : 05 Mar 2003
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ShortTermSqueezeTweak.hpp"
#include "edginc/Maths.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE

ShortTermSqueezeTweak::Shift::~Shift(){} // empty
ShortTermSqueezeTweak::RestorableShift::~RestorableShift(){} // empty

/** Sens Control for SHORT_TERM_SQUEEZE_TWEAK       */
const string ShortTermSqueezeTweak::NAME = "SHORT_TERM_SQUEEZE_TWEAK";
const double ShortTermSqueezeTweak::DEFAULT_SHIFT = 0.01;

/** constructor with explicit shift size */
ShortTermSqueezeTweak::ShortTermSqueezeTweak(double shiftSize):
    ScalarShift(TYPE, NAME, shiftSize), useShiftSizeSign(false) {
    validatePop2Object();
}
ShortTermSqueezeTweak::ShortTermSqueezeTweak(double shiftSize, bool useShiftSizeSign):
    ScalarShift(TYPE, NAME, shiftSize), useShiftSizeSign(useShiftSizeSign) {
    validatePop2Object();
}

/** for reflection */
ShortTermSqueezeTweak::ShortTermSqueezeTweak(): 
    ScalarShift(TYPE, NAME, DEFAULT_SHIFT), useShiftSizeSign(false){
}

/** Validation */
void ShortTermSqueezeTweak::validatePop2Object(){
    static const string method("ShortTermSqueezeTweak::validatePop2Object");
    double shiftSize = getShiftSize();
    if ( (shiftSize > 1.0) || (shiftSize < -1.0) ){
        throw ModelException( method, "ShortTermSqueezeTweak "+Format::toString(shiftSize)+
                              "must be between -1.0 and 1.0.");
    }
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double ShortTermSqueezeTweak::divisor() const{
    static const char routine[] = "ShortTermSqueezeTweak::divisor";
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
                         " in the SHORT_TERM_SQUEEZE_TWEAK calculation");
                throw e;
            }
        }
        // then just scale by shift size 
        double shiftSize = getShiftSize();    
        if (Maths::isZero(shiftSize)) {
            throw ModelException(routine, "Shift size is zero");
        }
        // report a 0.1 move in correlation squeeze (in [-1,1]); shift towards 0.0
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
CClassConstSP ShortTermSqueezeTweak::shiftInterface() const{
    return Shift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP ShortTermSqueezeTweak::restorableShiftInterface() const{
    return RestorableShift::TYPE;
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity.  The object must implement the
    VegaParallel.Shift interface */
bool ShortTermSqueezeTweak::nameMatches(const OutputName&         name,
                        IObjectConstSP            obj){
    // cast obj to ShortTermSqueezeTweak::Shift and then invoke name method
    const Shift& ShortTermSqueezeTweakObj = dynamic_cast<const Shift&>(*obj);
    return name.equals(ShortTermSqueezeTweakObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to
    this sensitivity to the supplied list */
void ShortTermSqueezeTweak::appendName(OutputNameArray&          namesList,
                       IObjectConstSP            obj){
    // cast obj to ShortTermSqueezeTweak::Shift and then invoke name method
    const Shift& ShortTermSqueezeTweakObj = dynamic_cast<const Shift&>(*obj);
    OutputNameSP outputName(new OutputName(ShortTermSqueezeTweakObj.sensName(this)));
    namesList.push_back(outputName);
}

bool ShortTermSqueezeTweak::shift(IObjectSP obj) {
    // cast obj to ShortTermSqueezeTweak::Shift and then invoke shift method
    Shift& ShortTermSqueezeTweakObj = dynamic_cast<Shift&>(*obj);
    return ShortTermSqueezeTweakObj.sensShift(this);
}

    
void ShortTermSqueezeTweak::restore(IObjectSP obj) {
    // cast obj to ShortTermSqueezeTweak::Shift and then invoke restore method
    RestorableShift& ShortTermSqueezeTweakObj = dynamic_cast<RestorableShift&>(*obj);
    ShortTermSqueezeTweakObj.sensRestore(this);
}

/** Returns true if implementations of Shift should override the sign of
    the shift size when making the tweak */
bool ShortTermSqueezeTweak::overrideShiftSizeSign() const{
    return !useShiftSizeSign;
}

class ShortTermSqueezeTweakHelper{
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new ShortTermSqueezeTweak(ShortTermSqueezeTweak::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new ShortTermSqueezeTweak(shiftSize);
        }
    };
    
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ShortTermSqueezeTweak, clazz);
        SUPERCLASS(ScalarShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultShortTermSqueezeTweak);
        FIELD_NO_DESC(useShiftSizeSign);
        FIELD_MAKE_TRANSIENT(useShiftSizeSign);
        // register how to build our sensitivity
        SensitivityFactory::addSens(ShortTermSqueezeTweak::NAME, 
                                    new Factory(), 
                                    new ShortTermSqueezeTweak(ShortTermSqueezeTweak::DEFAULT_SHIFT, false),
                                    ShortTermSqueezeTweak::Shift::TYPE);
    }

    static IObject* defaultShortTermSqueezeTweak(){
        return new ShortTermSqueezeTweak();
    }
};

CClassConstSP const ShortTermSqueezeTweak::TYPE = CClass::registerClassLoadMethod(
    "ShortTermSqueezeTweak", typeid(ShortTermSqueezeTweak), ShortTermSqueezeTweakHelper::load);

CClassConstSP const ShortTermSqueezeTweak::Shift::TYPE =
CClass::registerInterfaceLoadMethod(
    "ShortTermSqueezeTweak::Shift", typeid(ShortTermSqueezeTweak::Shift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(ShortTermSqueezeTweak::RestorableShift, clazz);
    EXTENDS(ShortTermSqueezeTweak::Shift);
}
    
CClassConstSP const ShortTermSqueezeTweak::RestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "ShortTermSqueezeTweak::RestorableShift", typeid(ShortTermSqueezeTweak::RestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE






