//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ContangoRhoPointwise.cpp
//
//   Description : Sensitivity to a pointwise shift in the contango curve
//
//   Author      : Andrew McCleery
//
//   Date        : 15 December 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ContangoRhoPointwise.hpp"
#include "edginc/Maths.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE
ContangoRhoPointwise::IShift::~IShift(){} // empty
ContangoRhoPointwise::IRestorableShift::~IRestorableShift(){} // empty

const double ONE_BASIS_POINT = 0.0001;

/** Sens Control for rho pointwise */
const string ContangoRhoPointwise::NAME = "CONTANGO_RHO_POINTWISE";
const double ContangoRhoPointwise::DEFAULT_SHIFT = ONE_BASIS_POINT;

/** constructor with explicit shift size */
ContangoRhoPointwise::ContangoRhoPointwise(double     shiftSize):
    VectorShift(TYPE, NAME, shiftSize){
}

/** for reflection */
ContangoRhoPointwise::ContangoRhoPointwise(): 
    VectorShift(TYPE, NAME){
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double ContangoRhoPointwise::divisor() const{
    static const string method = "ContangoRhoPointwise::divisor";
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
CClassConstSP ContangoRhoPointwise::shiftInterface() const{
    return IShift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP ContangoRhoPointwise::restorableShiftInterface() const{
    return IRestorableShift::TYPE;
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement the
    Shift interface */
bool ContangoRhoPointwise::nameMatches(const OutputName&         name,
                               IObjectConstSP            obj){
    // cast obj to ContangoRhoPointwise::Shift and then invoke name method
    const IShift& rhoPointwiseObj =
        dynamic_cast<const IShift&>(*obj);
    return name.equals(rhoPointwiseObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement the
    Shift interface */
void ContangoRhoPointwise::appendName(OutputNameArray&          namesList,
                              IObjectConstSP            obj){
    // cast obj to ContangoRhoPointwise::Shift and then invoke name method
    const IShift& rhoPointwiseObj =
        dynamic_cast<const IShift&>(*obj);
    OutputNameSP outputName(new OutputName(rhoPointwiseObj.sensName(this)));
    namesList.push_back(outputName);
}

bool ContangoRhoPointwise::shift(IObjectSP obj) {
    // cast obj to ContangoRhoPointwise::Shift and then invoke shift method
    IShift& rhoPointwiseObj =
        dynamic_cast<IShift&>(*obj);
    return rhoPointwiseObj.sensShift(this);
}

    
void ContangoRhoPointwise::restore(IObjectSP obj) {
    // cast obj to ContangoRhoPointwise::Shift and then invoke restore method
    IRestorableShift& rhoObj =
        dynamic_cast<IRestorableShift&>(*obj);
    rhoObj.sensRestore(this);
}

/** The supplied object is queried for the expiries array needed
    for doing rho pointwise and this array is returned. The supplied
    object must implement the ContangoRhoPointwise.Shift interface */
IObjectConstSP ContangoRhoPointwise::qualifier(IObjectConstSP obj){
    // cast obj to ContangoRhoPointwise::Shift and then invoke qualifier method
    const IShift& rhoObj =
        dynamic_cast<const IShift&>(*obj);
    return rhoObj.sensExpiries(this);
}

 
class ContangoRhoPointwiseHelper{
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new ContangoRhoPointwise(ContangoRhoPointwise::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new ContangoRhoPointwise(shiftSize);
        }
    };
        
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ContangoRhoPointwise, clazz);
        SUPERCLASS(VectorShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultContangoRhoPointwise);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(ContangoRhoPointwise::NAME, 
                                    new Factory(), 
                                    new ContangoRhoPointwise(ContangoRhoPointwise::DEFAULT_SHIFT),
                                    ContangoRhoPointwise::IShift::TYPE);
    }

    static IObject* defaultContangoRhoPointwise(){
        return new ContangoRhoPointwise();
    }
};

CClassConstSP const ContangoRhoPointwise::TYPE = CClass::registerClassLoadMethod(
    "ContangoRhoPointwise", typeid(ContangoRhoPointwise), ContangoRhoPointwiseHelper::load);

CClassConstSP const ContangoRhoPointwise::IShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "ContangoRhoPointwise::IShift", typeid(ContangoRhoPointwise::IShift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(ContangoRhoPointwise::IRestorableShift, clazz);
    EXTENDS(ContangoRhoPointwise::IShift);
}
    
CClassConstSP const ContangoRhoPointwise::IRestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "ContangoRhoPointwise::IRestorableShift",
    typeid(ContangoRhoPointwise::IRestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE
