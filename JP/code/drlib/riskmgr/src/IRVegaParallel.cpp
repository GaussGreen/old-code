//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IRVegaParallel.cpp
//
//   Description : Controls calculation of parallel IR vega
//
//   Author      : Andrew J Swain
//
//   Date        : 22 February 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/IRVegaParallel.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE
IRVegaParallel::Shift::~Shift(){} // empty
IRVegaParallel::RestorableShift::~RestorableShift(){} // empty

/** Sens Control for vega parallel */
const string IRVegaParallel::NAME = "IRVEGA_PARALLEL";
const double IRVegaParallel::DEFAULT_SHIFT = 0.001;

/** constructor with explicit shift size */
IRVegaParallel::IRVegaParallel(double shiftSize):
    ScalarShift(TYPE, NAME, shiftSize){
}

/** for reflection */
IRVegaParallel::IRVegaParallel(): 
    ScalarShift(TYPE, NAME){
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double IRVegaParallel::divisor() const{
    // add check for 0 shift size - somewhere 
    return 100.0 * getShiftSize();
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP IRVegaParallel::shiftInterface() const{
    return Shift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP IRVegaParallel::restorableShiftInterface() const{
    return RestorableShift::TYPE;
}

/** Returns true if the supplied object matches the supplied name
    for this sensitivity.  The object must implement the
    IRVegaParallel.Shift interface */
bool IRVegaParallel::nameMatches(const OutputName&         name,
                               IObjectConstSP          obj){
    // cast obj to IRVegaParallel::Shift and then invoke name method
    const Shift& vegaObj =
        dynamic_cast<const Shift&>(*obj);
    return name.equals(vegaObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to
    this sensitivity to the supplied list */
void IRVegaParallel::appendName(OutputNameArray&          namesList,
                              IObjectConstSP          obj){
    // cast obj to IRVegaParallel::Shift and then invoke name method
    const Shift& vegaObj =
        dynamic_cast<const Shift&>(*obj);
    OutputNameSP outputName(new OutputName(vegaObj.sensName(this)));
    namesList.push_back(outputName);
}

bool IRVegaParallel::shift(IObjectSP obj) {
    // cast obj to IRVegaParallel::Shift and then invoke shift method
    Shift& vegaObj =
        dynamic_cast<Shift&>(*obj);
    return vegaObj.sensShift(this);
}

void IRVegaParallel::restore(IObjectSP obj) {
    // cast obj to IRVegaParallel::Shift and then invoke restore method
    RestorableShift& vegaObj =
        dynamic_cast<RestorableShift&>(*obj);
    vegaObj.sensRestore(this);
}

class IRVegaParallelHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new IRVegaParallel(IRVegaParallel::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new IRVegaParallel(shiftSize);
        }
    };
        
public:
    /** Invoked when this class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(IRVegaParallel, clazz);
        SUPERCLASS(ScalarShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultIRVegaParallel);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(IRVegaParallel::NAME, 
                                    new Factory(), 
                                    new IRVegaParallel(IRVegaParallel::DEFAULT_SHIFT),
                                    IRVegaParallel::Shift::TYPE);
    }

    static IObject* defaultIRVegaParallel(){
        return new IRVegaParallel();
    }
};

CClassConstSP const IRVegaParallel::TYPE = CClass::registerClassLoadMethod(
    "IRVegaParallel", typeid(IRVegaParallel), IRVegaParallelHelper::load);

CClassConstSP const IRVegaParallel::Shift::TYPE =
CClass::registerInterfaceLoadMethod(
    "IRVegaParallel::Shift", typeid(IRVegaParallel::Shift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(IRVegaParallel::RestorableShift, clazz);
    EXTENDS(IRVegaParallel::Shift);
}
    
CClassConstSP const IRVegaParallel::RestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "IRVegaParallel::RestorableShift", typeid(IRVegaParallel::RestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE
