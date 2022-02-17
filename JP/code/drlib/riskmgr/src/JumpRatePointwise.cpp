//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : JumpRatePointwise.cpp
//
//   Description : for JumpRate pointwise
//
//   Date        : 25 Feb 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/JumpRatePointwise.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE

JumpRatePointwise::IShift::~IShift(){} // empty
JumpRatePointwise::IRestorableShift::~IRestorableShift(){} // empty

/** Sens Control for jump rate pointwise */
const string JumpRatePointwise::NAME = "JUMP_RATE_POINTWISE";
const double JumpRatePointwise::DEFAULT_SHIFT = 0.1;

/** constructor with explicit shift size */
JumpRatePointwise::JumpRatePointwise(double     shiftSize):
    VectorShift(TYPE, NAME, shiftSize){
}

/** constructor with explicit shift size and name override */ 
JumpRatePointwise::JumpRatePointwise(const string& overrideName,
                             double     shiftSize):
    VectorShift(TYPE, overrideName, shiftSize){
}

/** for reflection */
JumpRatePointwise::JumpRatePointwise(): 
    VectorShift(TYPE, NAME){
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double JumpRatePointwise::divisor() const{
    // add check for 0 shift size - somewhere 
    return 100.0 * getShiftSize();
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP JumpRatePointwise::shiftInterface() const{
    return IShift::TYPE;
}
 
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool JumpRatePointwise::nameMatches(const OutputName&         name,
                        IObjectConstSP          obj){
    // cast obj to JumpRatePointwise::Shift and then invoke name method
    const IShift& jumpRateObj =
        dynamic_cast<const IShift&>(*obj);
    return name.equals(jumpRateObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's Shift interface */
void JumpRatePointwise::appendName(OutputNameArray&          namesList,
                               IObjectConstSP          obj){
    // cast obj to JumpRatePointwise::Shift and then invoke name method
    const IShift& jumpRateObj =
        dynamic_cast<const IShift&>(*obj);
    OutputNameSP outputName(new OutputName(jumpRateObj.sensName(this)));
    namesList.push_back(outputName);
}

bool JumpRatePointwise::shift(IObjectSP obj) {
    // cast obj to JumpRatePointwise::Shift and then invoke shift method
    IShift& jumpRateObj =
        dynamic_cast<IShift&>(*obj);
    return jumpRateObj.sensShift(this);
}

/** The supplied object is queried for the expiries array needed
    for doing jump rate pointwise and this array is returned. The supplied
    object must implement the JumpRatePointwise.Shift interface */
IObjectConstSP JumpRatePointwise::qualifier(IObjectConstSP obj){
    // cast obj to JumpRatePointwise::Shift and then invoke qualifier method
    const IShift& jumpRateObj =
        dynamic_cast<const IShift&>(*obj);
    return jumpRateObj.sensExpiries(this);
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP JumpRatePointwise::restorableShiftInterface() const{
    return IRestorableShift::TYPE;
}

void JumpRatePointwise::restore(IObjectSP obj) {
    // cast obj to JumpRatePointwise::Shift and then invoke restore method
    IRestorableShift& jumpRateObj =
        dynamic_cast<IRestorableShift&>(*obj);
    jumpRateObj.sensRestore(this);
}

class JumpRatePointwiseHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new JumpRatePointwise(JumpRatePointwise::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new JumpRatePointwise(shiftSize);
        }
    };
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(JumpRatePointwise, clazz);
        SUPERCLASS(VectorShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultJumpRatePointwise);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(JumpRatePointwise::NAME, 
                                    new Factory(), 
                                    new JumpRatePointwise(JumpRatePointwise::DEFAULT_SHIFT),
                                    JumpRatePointwise::IShift::TYPE);
    }

    static IObject* defaultJumpRatePointwise(){
        return new JumpRatePointwise();
    }
};

CClassConstSP const JumpRatePointwise::TYPE = CClass::registerClassLoadMethod(
    "JumpRatePointwise", typeid(JumpRatePointwise), JumpRatePointwiseHelper::load);

CClassConstSP const JumpRatePointwise::IShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "JumpRatePointwise::IShift", typeid(JumpRatePointwise::IShift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(JumpRatePointwise::IRestorableShift, clazz);
    EXTENDS(JumpRatePointwise::IShift);
}
    
CClassConstSP const JumpRatePointwise::IRestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "JumpRatePointwise::IRestorableShift",
    typeid(JumpRatePointwise::IRestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE
