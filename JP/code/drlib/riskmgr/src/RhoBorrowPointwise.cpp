//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RhoBorrowPointwise.cpp
//
//   Description : Controls calculation of Rho pointwise
//
//   Author      : Stephen Hope
//
//   Date        : 9 March 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/RhoBorrowPointwise.hpp"
#include "edginc/Maths.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE
RhoBorrowPointwise::IShift::~IShift(){} // empty
RhoBorrowPointwise::IRestorableShift::~IRestorableShift(){} // empty

const double ONE_BASIS_POINT = 0.0001;

/** Sens Control for rho pointwise */
const string RhoBorrowPointwise::NAME = "RHO_BORROW_POINTWISE";
const double RhoBorrowPointwise::DEFAULT_SHIFT = ONE_BASIS_POINT;

/** constructor with explicit shift size */
RhoBorrowPointwise::RhoBorrowPointwise(double     shiftSize):
    VectorShift(TYPE, NAME, shiftSize){
}

/** for reflection */
RhoBorrowPointwise::RhoBorrowPointwise(): 
    VectorShift(TYPE, NAME){
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double RhoBorrowPointwise::divisor() const{
    static const string method = "RhoBorrowPointwise::divisor";
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
CClassConstSP RhoBorrowPointwise::shiftInterface() const{
    return IShift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP RhoBorrowPointwise::restorableShiftInterface() const{
    return IRestorableShift::TYPE;
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity.  The object must implement the
    VegaPointwise.Shift interface */
bool RhoBorrowPointwise::nameMatches(const OutputName&         name,
                        IObjectConstSP          obj){
    // cast obj to RhoBorrowPointwise::Shift and then invoke name method
    const IShift& rhoBorrowObj =
        dynamic_cast<const IShift&>(*obj);
    return name.equals(rhoBorrowObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to
    this sensitivity to the supplied list */
void RhoBorrowPointwise::appendName(OutputNameArray&          namesList,
                       IObjectConstSP          obj){
    // cast obj to RhoBorrowPointwise::Shift and then invoke name method
    const IShift& rhoBorrowObj =
        dynamic_cast<const IShift&>(*obj);
    OutputNameSP outputName(new OutputName(rhoBorrowObj.sensName(this)));
    namesList.push_back(outputName);
}

bool RhoBorrowPointwise::shift(IObjectSP obj) {
    // cast obj to RhoBorrowPointwise::Shift and then invoke shift method
    IShift& rhoObj =
        dynamic_cast<IShift&>(*obj);
    return rhoObj.sensShift(this);
}

bool RhoBorrowPointwise::restorableShift(IObjectConstSP& obj) const{
    return IRestorableShift::TYPE->isInstance(obj);
}
    
void RhoBorrowPointwise::restore(IObjectSP obj) {
    // cast obj to RhoBorrowPointwise::Shift and then invoke restore method
    IRestorableShift& rhoObj =
        dynamic_cast<IRestorableShift&>(*obj);
    rhoObj.sensRestore(this);
}

/** The supplied object is queried for the expiries array needed
    for doing rho pointwise and this array is returned. The supplied
    object must implement the RhoBorrowPointwise.Shift interface */
IObjectConstSP RhoBorrowPointwise::qualifier(IObjectConstSP obj){
    // cast obj to RhoBorrowPointwise::Shift and then invoke qualifier method
    const IShift& rhoObj =
        dynamic_cast<const IShift&>(*obj);
    return rhoObj.sensExpiries(this);
}

 
class RhoBorrowPointwiseHelper{
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new RhoBorrowPointwise(RhoBorrowPointwise::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new RhoBorrowPointwise(shiftSize);
        }
    };
        
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(RhoBorrowPointwise, clazz);
        SUPERCLASS(VectorShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultRhoBorrowPointwise);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(RhoBorrowPointwise::NAME, 
                                    new Factory(), 
                                    new RhoBorrowPointwise(RhoBorrowPointwise::DEFAULT_SHIFT),
                                    RhoBorrowPointwise::IShift::TYPE);
    }

    static IObject* defaultRhoBorrowPointwise(){
        return new RhoBorrowPointwise();
    }
};

CClassConstSP const RhoBorrowPointwise::TYPE = CClass::registerClassLoadMethod(
    "RhoBorrowPointwise", typeid(RhoBorrowPointwise), RhoBorrowPointwiseHelper::load);

CClassConstSP const RhoBorrowPointwise::IShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "RhoBorrowPointwise::IShift", typeid(RhoBorrowPointwise::IShift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(RhoBorrowPointwise::IRestorableShift, clazz);
    EXTENDS(RhoBorrowPointwise::IShift);
}
    
CClassConstSP const RhoBorrowPointwise::IRestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "RhoBorrowPointwise::IRestorableShift",
    typeid(RhoBorrowPointwise::IRestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE
