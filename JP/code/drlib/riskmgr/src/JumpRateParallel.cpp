//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : JumpRateParallel.cpp
//
//   Description : for JumpRate parallel tweak
//
//   Date        : 25 Feb 2002
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/JumpRateParallel.hpp"

DRLIB_BEGIN_NAMESPACE

JumpRateParallel::Shift::~Shift(){} // empty
JumpRateParallel::RestorableShift::~RestorableShift(){} // empty

/** Sens Control for jump rate parallel */
const string JumpRateParallel::NAME = "JUMP_RATE_PARALLEL";
const double JumpRateParallel::DEFAULT_SHIFT = 0.1;

/** constructor with explicit shift size */
JumpRateParallel::JumpRateParallel(double     shiftSize):
    ScalarShift(TYPE, NAME, shiftSize){
}

/** for reflection */
JumpRateParallel::JumpRateParallel(): 
    ScalarShift(TYPE, NAME){
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double JumpRateParallel::divisor() const{
    // add check for 0 shift size - somewhere 
    return 100.0 * getShiftSize();
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP JumpRateParallel::shiftInterface() const{
    return Shift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP JumpRateParallel::restorableShiftInterface() const{
    return RestorableShift::TYPE;
}

/** Returns true if the supplied object matches the supplied name
    for this sensitivity.  The object must implement the
    JumpRateParallel.Shift interface */
bool JumpRateParallel::nameMatches(const OutputName&         name,
                               IObjectConstSP          obj){
    // cast obj to JumpRateParallel::Shift and then invoke name method
    const Shift& jumpRateObj =
        dynamic_cast<const Shift&>(*obj);
    return name.equals(jumpRateObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to
    this sensitivity to the supplied list */
void JumpRateParallel::appendName(OutputNameArray&          namesList,
                              IObjectConstSP          obj){
    // cast obj to JumpRateParallel::Shift and then invoke name method
    const Shift& jumpRateObj =
        dynamic_cast<const Shift&>(*obj);
    OutputNameSP outputName(new OutputName(jumpRateObj.sensName(this)));
    namesList.push_back(outputName);
}

bool JumpRateParallel::shift(IObjectSP obj) {
    // cast obj to JumpRateParallel::Shift and then invoke shift method
    Shift& jumpRateObj =
        dynamic_cast<Shift&>(*obj);
    return jumpRateObj.sensShift(this);
}

void JumpRateParallel::restore(IObjectSP obj) {
    // cast obj to JumpRateParallel::Shift and then invoke restore method
    RestorableShift& jumpRateObj =
        dynamic_cast<RestorableShift&>(*obj);
    jumpRateObj.sensRestore(this);
}

class JumpRateParallelHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new JumpRateParallel(JumpRateParallel::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new JumpRateParallel(shiftSize);
        }
    };
        
public:
    /** Invoked when this class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(JumpRateParallel, clazz);
        SUPERCLASS(ScalarShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultJumpRateParallel);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(JumpRateParallel::NAME, 
                                    new Factory(), 
                                    new JumpRateParallel(JumpRateParallel::DEFAULT_SHIFT),
                                    JumpRateParallel::Shift::TYPE);
    }

    static IObject* defaultJumpRateParallel(){
        return new JumpRateParallel();
    }
};

CClassConstSP const JumpRateParallel::TYPE = CClass::registerClassLoadMethod(
    "JumpRateParallel", typeid(JumpRateParallel), JumpRateParallelHelper::load);

CClassConstSP const JumpRateParallel::Shift::TYPE =
CClass::registerInterfaceLoadMethod(
    "JumpRateParallel::Shift", typeid(JumpRateParallel::Shift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(JumpRateParallel::RestorableShift, clazz);
    EXTENDS(JumpRateParallel::Shift);
}
    
CClassConstSP const JumpRateParallel::RestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "JumpRateParallel::RestorableShift", typeid(JumpRateParallel::RestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE
