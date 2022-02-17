//----------------------------------------------------------------------------
//
//   Group       : Energy Exotics Derivatives Research
//
//   Filename    : ForwardParallel.cpp
//
//   Description : Sensitivity to a parallel shift in the forward curve
//
//   Author      : Simon Creeger
//
//   Date        : 12 September 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ForwardParallel.hpp"
#include "edginc/Maths.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE
ForwardParallel::Shift::~Shift(){} // empty
ForwardParallel::RestorableShift::~RestorableShift(){} // empty

/** Sens Control for forward parallel */
const string ForwardParallel::NAME = "FORWARD_PARALLEL";
const double ForwardParallel::DEFAULT_SHIFT = 1.0;
 
/** constructor with explicit shift size */
ForwardParallel::ForwardParallel(double shiftSize):
    ScalarShift(TYPE, NAME, shiftSize){
}

/** constructor with explicit shift size */
ForwardParallel::ForwardParallel(double     shiftSize,
             IModel*    model,
             Control*   control):
    ScalarShift(TYPE, NAME, shiftSize) {
    this->algorithm = model;
    this->control   = control;
}

/** for reflection */
ForwardParallel::ForwardParallel(): 
    ScalarShift(TYPE, NAME){
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double ForwardParallel::divisor() const{
    static const string method = "ForwardParallel::divisor";
    double shiftSize;
    try{
        // just scale the shift size for a 1bp move
        shiftSize = getShiftSize();
        if (Maths::isZero(shiftSize)){
            throw ModelException(method, "Shift size is zero");
        }
        return shiftSize;
    } 
    catch (ModelException& e) {
        throw ModelException(&e, method);
    }
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP ForwardParallel::shiftInterface() const{
    return Shift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP ForwardParallel::restorableShiftInterface() const{
    return RestorableShift::TYPE;
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement the
    Shift interface */
bool ForwardParallel::nameMatches(const OutputName&         name,
                                      IObjectConstSP            obj){
    // cast obj to ForwardParallel::Shift and then invoke name method
    const Shift& forwardParallelObj =
        dynamic_cast<const Shift&>(*obj);
    return name.equals(forwardParallelObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement the
    Shift interface */
void ForwardParallel::appendName(OutputNameArray&          namesList,
                              IObjectConstSP            obj){
    // cast obj to ForwardParallel::Shift and then invoke name method
    const Shift& forwardParallelObj =
        dynamic_cast<const Shift&>(*obj);
    OutputNameSP outputName(new OutputName(forwardParallelObj.sensName(this)));
    namesList.push_back(outputName);
}

bool ForwardParallel::shift(IObjectSP obj) {
    // cast obj to ForwardParallel::Shift and then invoke shift method
    Shift& forwardParallelObj =
        dynamic_cast<Shift&>(*obj);
    return forwardParallelObj.sensShift(this);
}

    
void ForwardParallel::restore(IObjectSP obj) {
    // cast obj to ForwardParallel::Shift and then invoke restore method
    RestorableShift& forwardParallelObj =
        dynamic_cast<RestorableShift&>(*obj);
    forwardParallelObj.sensRestore(this);
}

class ForwardParallelHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new ForwardParallel(ForwardParallel::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new ForwardParallel(shiftSize);
        }
    };
        
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ForwardParallel, clazz);
        SUPERCLASS(ScalarShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultForwardParallel);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(ForwardParallel::NAME, 
                                    new Factory(), 
                                    new ForwardParallel(ForwardParallel::DEFAULT_SHIFT),
                                    ForwardParallel::Shift::TYPE);
    }

    static IObject* defaultForwardParallel(){
        return new ForwardParallel();
    }
};

CClassConstSP const ForwardParallel::TYPE = CClass::registerClassLoadMethod(
    "ForwardParallel", typeid(ForwardParallel), ForwardParallelHelper::load);

CClassConstSP const ForwardParallel::Shift::TYPE =
CClass::registerInterfaceLoadMethod(
    "ForwardParallel::Shift", typeid(ForwardParallel::Shift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(ForwardParallel::RestorableShift, clazz);
    EXTENDS(ForwardParallel::Shift);
}
    
CClassConstSP const ForwardParallel::RestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "ForwardParallel::RestorableShift", typeid(ForwardParallel::RestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE
