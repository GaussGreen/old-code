//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ContangoRhoParallel.cpp
//
//   Description : Sensitivity to a parallel shift in the contango curve
//
//   Author      : Andrew McCleery
//
//   Date        : 15 December 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ContangoRhoParallel.hpp"
#include "edginc/Maths.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE
ContangoRhoParallel::Shift::~Shift(){} // empty
ContangoRhoParallel::RestorableShift::~RestorableShift(){} // empty

const double ONE_BASIS_POINT = 0.0001;

/** Sens Control for rho parallel */
const string ContangoRhoParallel::NAME = "CONTANGO_RHO_PARALLEL";
const double ContangoRhoParallel::DEFAULT_SHIFT = ONE_BASIS_POINT;

/** constructor with explicit shift size */
ContangoRhoParallel::ContangoRhoParallel(double shiftSize):
    ScalarShift(TYPE, NAME, shiftSize){
}

/** constructor with explicit shift size */
ContangoRhoParallel::ContangoRhoParallel(double     shiftSize,
             IModel*    model,
             Control*   control):
    ScalarShift(TYPE, NAME, shiftSize) {
    this->algorithm = model;
    this->control   = control;
}

/** for reflection */
ContangoRhoParallel::ContangoRhoParallel(): 
    ScalarShift(TYPE, NAME){
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double ContangoRhoParallel::divisor() const{
    static const string method = "ContangoRhoParallel::divisor";
    double shiftSize;
    try{
        // just scale the shift size for a 1bp move
        shiftSize = getShiftSize();
        if (Maths::isZero(shiftSize)){
            throw ModelException(method, "Shift size is zero");
        }
        return (shiftSize/ONE_BASIS_POINT);
    } 
    catch (ModelException& e) {
        throw ModelException(&e, method);
    }
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP ContangoRhoParallel::shiftInterface() const{
    return Shift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP ContangoRhoParallel::restorableShiftInterface() const{
    return RestorableShift::TYPE;
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement the
    Shift interface */
bool ContangoRhoParallel::nameMatches(const OutputName&         name,
                                      IObjectConstSP            obj){
    // cast obj to ContangoRhoParallel::Shift and then invoke name method
    const Shift& rhoParallelObj =
        dynamic_cast<const Shift&>(*obj);
    return name.equals(rhoParallelObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement the
    Shift interface */
void ContangoRhoParallel::appendName(OutputNameArray&          namesList,
                              IObjectConstSP            obj){
    // cast obj to ContangoRhoParallel::Shift and then invoke name method
    const Shift& rhoParallelObj =
        dynamic_cast<const Shift&>(*obj);
    OutputNameSP outputName(new OutputName(rhoParallelObj.sensName(this)));
    namesList.push_back(outputName);
}

bool ContangoRhoParallel::shift(IObjectSP obj) {
    // cast obj to ContangoRhoParallel::Shift and then invoke shift method
    Shift& rhoParallelObj =
        dynamic_cast<Shift&>(*obj);
    return rhoParallelObj.sensShift(this);
}

    
void ContangoRhoParallel::restore(IObjectSP obj) {
    // cast obj to ContangoRhoParallel::Shift and then invoke restore method
    RestorableShift& rhoParallelObj =
        dynamic_cast<RestorableShift&>(*obj);
    rhoParallelObj.sensRestore(this);
}

class ContangoRhoParallelHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new ContangoRhoParallel(ContangoRhoParallel::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new ContangoRhoParallel(shiftSize);
        }
    };
        
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ContangoRhoParallel, clazz);
        SUPERCLASS(ScalarShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultContangoRhoParallel);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(ContangoRhoParallel::NAME, 
                                    new Factory(), 
                                    new ContangoRhoParallel(ContangoRhoParallel::DEFAULT_SHIFT),
                                    ContangoRhoParallel::Shift::TYPE);
    }

    static IObject* defaultContangoRhoParallel(){
        return new ContangoRhoParallel();
    }
};

CClassConstSP const ContangoRhoParallel::TYPE = CClass::registerClassLoadMethod(
    "ContangoRhoParallel", typeid(ContangoRhoParallel), ContangoRhoParallelHelper::load);

CClassConstSP const ContangoRhoParallel::Shift::TYPE =
CClass::registerInterfaceLoadMethod(
    "ContangoRhoParallel::Shift", typeid(ContangoRhoParallel::Shift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(ContangoRhoParallel::RestorableShift, clazz);
    EXTENDS(ContangoRhoParallel::Shift);
}
    
CClassConstSP const ContangoRhoParallel::RestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "ContangoRhoParallel::RestorableShift", typeid(ContangoRhoParallel::RestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE
