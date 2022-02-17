//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RhoBorrowParallel.cpp
//
//   Description : Rho borrow parallel sensitivity
//
//   Author      : Stephen Hope
//
//   Date        : 26 February 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/RhoBorrowParallel.hpp"
#include "edginc/Maths.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE

RhoBorrowParallel::Shift::~Shift(){} // empty
RhoBorrowParallel::RestorableShift::~RestorableShift(){} // empty

const double ONE_BASIS_POINT = 0.0001;

/** Sens Control for rho parallel */
const string RhoBorrowParallel::NAME = "RHO_BORROW_PARALLEL";
const double RhoBorrowParallel::DEFAULT_SHIFT = ONE_BASIS_POINT;

/** constructor with explicit shift size */
RhoBorrowParallel::RhoBorrowParallel(double shiftSize):
    ScalarShift(TYPE, NAME, shiftSize){
}

/** for reflection */
RhoBorrowParallel::RhoBorrowParallel(): 
    ScalarShift(TYPE, NAME){
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double RhoBorrowParallel::divisor() const{
    static const string method = "RhoBorrowParallel::divisor";
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
CClassConstSP RhoBorrowParallel::shiftInterface() const{
    return Shift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP RhoBorrowParallel::restorableShiftInterface() const{
    return RestorableShift::TYPE;
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity.  The object must implement the
    VegaParallel.Shift interface */
bool RhoBorrowParallel::nameMatches(const OutputName&         name,
                                    IObjectConstSP            obj){
    // cast obj to RhoBorrowParallel::Shift and then invoke name method
    const Shift& rhoBorrowObj = dynamic_cast<const Shift&>(*obj);
    return name.equals(rhoBorrowObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to
    this sensitivity to the supplied list */
void RhoBorrowParallel::appendName(OutputNameArray&          namesList,
                                   IObjectConstSP            obj){
    // cast obj to RhoBorrowParallel::Shift and then invoke name method
    const Shift& rhoBorrowObj = dynamic_cast<const Shift&>(*obj);
    OutputNameSP outputName(new OutputName(rhoBorrowObj.sensName(this)));
    namesList.push_back(outputName);
}

bool RhoBorrowParallel::shift(IObjectSP obj) {
    // cast obj to RhoBorrowParallel::Shift and then invoke shift method
    Shift& rhoBorrowObj = dynamic_cast<Shift&>(*obj);
    return rhoBorrowObj.sensShift(this);
}

    
void RhoBorrowParallel::restore(IObjectSP obj) {
    // cast obj to RhoBorrowParallel::Shift and then invoke restore method
    RestorableShift& rhoBorrowObj = dynamic_cast<RestorableShift&>(*obj);
    rhoBorrowObj.sensRestore(this);
}

class RhoBorrowParallelHelper{
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new RhoBorrowParallel(RhoBorrowParallel::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new RhoBorrowParallel(shiftSize);
        }
    };
        
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(RhoBorrowParallel, clazz);
        SUPERCLASS(ScalarShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultRhoBorrowParallel);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(RhoBorrowParallel::NAME, 
                                    new Factory(), 
                                    new RhoBorrowParallel(RhoBorrowParallel::DEFAULT_SHIFT),
                                    RhoBorrowParallel::Shift::TYPE);
    }

    static IObject* defaultRhoBorrowParallel(){
        return new RhoBorrowParallel();
    }
};

CClassConstSP const RhoBorrowParallel::TYPE = CClass::registerClassLoadMethod(
    "RhoBorrowParallel", typeid(RhoBorrowParallel), RhoBorrowParallelHelper::load);

CClassConstSP const RhoBorrowParallel::Shift::TYPE =
CClass::registerInterfaceLoadMethod(
    "RhoBorrowParallel::Shift", typeid(RhoBorrowParallel::Shift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(RhoBorrowParallel::RestorableShift, clazz);
    EXTENDS(RhoBorrowParallel::Shift);
}
    
CClassConstSP const RhoBorrowParallel::RestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "RhoBorrowParallel::RestorableShift", typeid(RhoBorrowParallel::RestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE
