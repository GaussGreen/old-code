//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MuParallel.cpp
//
//   Description : MU_PARALLEL sensitivity
//
//   Author      : Stephen Hope
//
//   Date        : 12 March 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/MuParallel.hpp"
#include "edginc/Maths.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE
MuParallel::Shift::~Shift(){} // empty
MuParallel::RestorableShift::~RestorableShift(){} // empty

/** Sens Control for MU_PARALLEL */
const string MuParallel::NAME = "MU_PARALLEL";
const double MuParallel::DEFAULT_SHIFT = 0.1;

/** constructor with explicit shift size */
MuParallel::MuParallel(double shiftSize):
    ScalarShift(TYPE, NAME, shiftSize){
}

/** for reflection */
MuParallel::MuParallel(): 
    ScalarShift(TYPE, NAME){
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double MuParallel::divisor() const{
    static const char routine[] = "MuParallel::divisor";
    double shiftSize;
    try{
        shiftSize = getShiftSize();
        if (Maths::isZero(shiftSize) ){
            throw ModelException(routine, "Shift size is zero");
        }
    } catch (ModelException& e){
        e.addMsg(routine);
        throw e;
    }
    return (shiftSize * 100.0);
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP MuParallel::shiftInterface() const{
    return Shift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP MuParallel::restorableShiftInterface() const{
    return RestorableShift::TYPE;
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this sensitivity's
    Shift interface */
bool MuParallel::nameMatches(const OutputName&         name,
                               IObjectConstSP          obj){
    // cast obj to MuParallel::Shift and then invoke name method
    const Shift& muParallelObj = 
        dynamic_cast<const Shift&>(*obj);
    return name.equals(muParallelObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    sensitivity's Shift interface */
void MuParallel::appendName(OutputNameArray&          namesList,
                              IObjectConstSP          obj){
    // cast obj to MuParallel::Shift and then invoke name method
    const Shift& muParallelObj = 
        dynamic_cast<const Shift&>(*obj);
    OutputNameSP outputName(new OutputName(muParallelObj.sensName(this)));
    namesList.push_back(outputName);
}

bool MuParallel::shift(IObjectSP obj) {
    // cast obj to MuParallel::Shift and then invoke shift method
    Shift& muParallelObj = 
        dynamic_cast<Shift&>(*obj);
    return muParallelObj.sensShift(this);
}

void MuParallel::restore(IObjectSP obj) {
    // cast obj to MuParallel::Shift and then invoke restore method
    RestorableShift& muParallelObj = 
        dynamic_cast<RestorableShift&>(*obj);
    muParallelObj.sensRestore(this);
}

class MuParallelHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new MuParallel(MuParallel::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new MuParallel(shiftSize);
        }
    };
        
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(MuParallel, clazz);
        SUPERCLASS(ScalarShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultMuParallel);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(MuParallel::NAME, 
                                    new Factory(), 
                                    new MuParallel(MuParallel::DEFAULT_SHIFT),
                                    MuParallel::Shift::TYPE);
    }

    static IObject* defaultMuParallel(){
        return new MuParallel();
    }
};

CClassConstSP const MuParallel::TYPE = CClass::registerClassLoadMethod(
    "MuParallel", typeid(MuParallel), MuParallelHelper::load);

CClassConstSP const MuParallel::Shift::TYPE =
CClass::registerInterfaceLoadMethod(
    "MuParallel::Shift", typeid(MuParallel::Shift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(MuParallel::RestorableShift, clazz);
    EXTENDS(MuParallel::Shift);
}
    
CClassConstSP const MuParallel::RestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "MuParallel::RestorableShift", typeid(MuParallel::RestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE
