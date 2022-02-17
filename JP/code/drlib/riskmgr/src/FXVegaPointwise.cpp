//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FXVegaPointwise.cpp
//
//   Description : Controls calculation of FX Vega pointwise
//
//   Author      : Andrew J Swain
//
//   Date        : 3 February 2003
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/FXVegaPointwise.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE
FXVegaPointwise::IShift::~IShift(){} // empty
FXVegaPointwise::IRestorableShift::~IRestorableShift(){} // empty

/** Sens Control for vega pointwise */
const string FXVegaPointwise::NAME = "FX_VEGA_POINTWISE";
const double FXVegaPointwise::DEFAULT_SHIFT = 0.001;

/** constructor with explicit shift size */
FXVegaPointwise::FXVegaPointwise(double     shiftSize):
    VectorShift(TYPE, NAME, shiftSize){
}

/** constructor with explicit shift size and explicit expiry to
    tweak. Also need to specify all the expiries. Building
    FXVegaPointwise this way is useful as it allows individal points to
    be tweaked. Note that the supplied expiries are overwritten by the
    vector shift calculate method */
FXVegaPointwise::FXVegaPointwise(double                     shiftSize,
                                 const ExpiryConstSP&       expiry,
                                 const ExpiryArrayConstSP&  allExpiries):
    VectorShift(TYPE, NAME, shiftSize){
    this->expiry = expiry; //expiry in parent
    this->cachedExpiries = allExpiries; //expiry in parent
}

/** constructor with explicit shift size and name override (override allows
    a VEGA_POINTWISE calculation to be stored under, eg, VEGA_MATRIX) */
FXVegaPointwise::FXVegaPointwise(const string& overrideName,
                                 double     shiftSize):
    VectorShift(TYPE, overrideName, shiftSize){
}

/** for reflection */
FXVegaPointwise::FXVegaPointwise(): 
    VectorShift(TYPE, NAME){
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double FXVegaPointwise::divisor() const{
    // add check for 0 shift size - somewhere 
    return 100.0 * getShiftSize();
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP FXVegaPointwise::shiftInterface() const{
    return IShift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP FXVegaPointwise::restorableShiftInterface() const{
    return IRestorableShift::TYPE;
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool FXVegaPointwise::nameMatches(const OutputName&         name,
                                  IObjectConstSP          obj){
    // cast obj to FXVegaPointwise::Shift and then invoke name method
    const IShift& vegaObj =
        dynamic_cast<const IShift&>(*obj);
    return name.equals(vegaObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's Shift interface */
void FXVegaPointwise::appendName(OutputNameArray&          namesList,
                                 IObjectConstSP          obj){
    // cast obj to FXVegaPointwise::Shift and then invoke name method
    const IShift& vegaObj = dynamic_cast<const IShift&>(*obj);
    OutputNameSP outputName(new OutputName(vegaObj.sensName(this)));
    namesList.push_back(outputName);
}

bool FXVegaPointwise::shift(IObjectSP obj) {
    // cast obj to FXVegaPointwise::Shift and then invoke shift method
    IShift& vegaObj = dynamic_cast<IShift&>(*obj);
    return vegaObj.sensShift(this);
}

void FXVegaPointwise::restore(IObjectSP obj) {
    // cast obj to FXVegaPointwise::Shift and then invoke restore method
    IRestorableShift& vegaObj = dynamic_cast<IRestorableShift&>(*obj);
    vegaObj.sensRestore(this);
}

/** The supplied object is queried for the expiries array needed
    for doing vega pointwise and this array is returned. The supplied
    object must implement the FXVegaPointwise.Shift interface */
IObjectConstSP FXVegaPointwise::qualifier(IObjectConstSP obj){
    // cast obj to FXVegaPointwise::Shift and then invoke qualifier method
    const IShift& vegaObj = dynamic_cast<const IShift&>(*obj);
    return vegaObj.sensExpiries(this);
}

class FXVegaPointwiseHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new FXVegaPointwise(FXVegaPointwise::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new FXVegaPointwise(shiftSize);
        }
    };
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(FXVegaPointwise, clazz);
        SUPERCLASS(VectorShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultFXVegaPointwise);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(FXVegaPointwise::NAME, 
                                    new Factory(), 
                                    new FXVegaPointwise(FXVegaPointwise::DEFAULT_SHIFT),
                                    FXVegaPointwise::IShift::TYPE);
    }

    static IObject* defaultFXVegaPointwise(){
        return new FXVegaPointwise();
    }
};

CClassConstSP const FXVegaPointwise::TYPE = CClass::registerClassLoadMethod(
    "FXVegaPointwise", typeid(FXVegaPointwise), FXVegaPointwiseHelper::load);

CClassConstSP const FXVegaPointwise::IShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "FXVegaPointwise::IShift", typeid(FXVegaPointwise::IShift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(FXVegaPointwise::IRestorableShift, clazz);
    EXTENDS(FXVegaPointwise::IShift);
}
    
CClassConstSP const FXVegaPointwise::IRestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "FXVegaPointwise::IRestorableShift",
    typeid(FXVegaPointwise::IRestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE
