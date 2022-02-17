//----------------------------------------------------------------------------
//
//   Group       : Energy Exotics Derivatives Research
//
//   Filename    : ForwardPointwise.cpp
//
//   Description : Sensitivity to a pointwise shift in the Forward curve
//
//   Author      : Simon Creeger
//
//   Date        : 12 September 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ForwardPointwise.hpp"
#include "edginc/Maths.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE
ForwardPointwise::IShift::~IShift(){} // empty
ForwardPointwise::IRestorableShift::~IRestorableShift(){} // empty

/** Sens Control for forward pointwise */
const string ForwardPointwise::NAME = "FORWARD_POINTWISE";
const double ForwardPointwise::DEFAULT_SHIFT = 1.0; // need to change this SAC

/** constructor with explicit shift size */
ForwardPointwise::ForwardPointwise(double     shiftSize):
    VectorShift(TYPE, NAME, shiftSize){
}

/** for reflection */
ForwardPointwise::ForwardPointwise(): 
    VectorShift(TYPE, NAME){
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double ForwardPointwise::divisor() const{
    static const string method = "ForwardPointwise::divisor";
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
CClassConstSP ForwardPointwise::shiftInterface() const{
    return IShift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP ForwardPointwise::restorableShiftInterface() const{
    return IRestorableShift::TYPE;
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement the
    Shift interface */
bool ForwardPointwise::nameMatches(const OutputName&         name,
                               IObjectConstSP            obj){
    // cast obj to ForwardPointwise::Shift and then invoke name method
    const IShift& forwardPointwiseObj =
        dynamic_cast<const IShift&>(*obj);
    return name.equals(forwardPointwiseObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement the
    Shift interface */
void ForwardPointwise::appendName(OutputNameArray&          namesList,
                              IObjectConstSP            obj){
    // cast obj to ForwardPointwise::Shift and then invoke name method
    const IShift& forwardPointwiseObj =
        dynamic_cast<const IShift&>(*obj);
    OutputNameSP outputName(new OutputName(forwardPointwiseObj.sensName(this)));
    namesList.push_back(outputName);
}

bool ForwardPointwise::shift(IObjectSP obj) {
    // cast obj to ForwardPointwise::Shift and then invoke shift method
    IShift& forwardPointwiseObj =
        dynamic_cast<IShift&>(*obj);
    return forwardPointwiseObj.sensShift(this);
}

    
void ForwardPointwise::restore(IObjectSP obj) {
    // cast obj to ForwardPointwise::Shift and then invoke restore method
    IRestorableShift& forwardObj =
        dynamic_cast<IRestorableShift&>(*obj);
    forwardObj.sensRestore(this);
}

/** The supplied object is queried for the expiries array needed
    for doing forward pointwise and this array is returned. The supplied
    object must implement the ForwardPointwise.Shift interface */
IObjectConstSP ForwardPointwise::qualifier(IObjectConstSP obj){
    // cast obj to ForwardPointwise::Shift and then invoke qualifier method
    const IShift& forwardObj =
        dynamic_cast<const IShift&>(*obj);
    return forwardObj.sensExpiries(this);
}

 
class ForwardPointwiseHelper{
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new ForwardPointwise(ForwardPointwise::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new ForwardPointwise(shiftSize);
        }
    };
        
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ForwardPointwise, clazz);
        SUPERCLASS(VectorShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultForwardPointwise);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(ForwardPointwise::NAME, 
                                    new Factory(), 
                                    new ForwardPointwise(ForwardPointwise::DEFAULT_SHIFT),
                                    ForwardPointwise::IShift::TYPE);
    }

    static IObject* defaultForwardPointwise(){
        return new ForwardPointwise();
    }
};

CClassConstSP const ForwardPointwise::TYPE = CClass::registerClassLoadMethod(
    "ForwardPointwise", typeid(ForwardPointwise), ForwardPointwiseHelper::load);

CClassConstSP const ForwardPointwise::IShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "ForwardPointwise::IShift", typeid(ForwardPointwise::IShift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(ForwardPointwise::IRestorableShift, clazz);
    EXTENDS(ForwardPointwise::IShift);
}
    
CClassConstSP const ForwardPointwise::IRestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "ForwardPointwise::IRestorableShift",
    typeid(ForwardPointwise::IRestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE
