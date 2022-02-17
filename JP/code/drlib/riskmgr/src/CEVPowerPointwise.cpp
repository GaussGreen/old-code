//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CEVPowerPointwise.cpp
//
//   Description : for CEVPower pointwise
//
//   Date        : 25 Feb 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CEVPowerPointwise.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE

// not sure why we have to have these - MSVC gets upset if they're not in the
// header
CEVPowerPointwise::IShift::~IShift(){} // empty
CEVPowerPointwise::IRestorableShift::~IRestorableShift(){} // empty

/** Sens Control for CEVPower pointwise */
const string CEVPowerPointwise::NAME = "CEVPOWER_POINTWISE";
const double CEVPowerPointwise::DEFAULT_SHIFT = 0.1;

/** constructor with explicit shift size */
CEVPowerPointwise::CEVPowerPointwise(double     shiftSize):
    VectorShift(TYPE, NAME, shiftSize){
}

/** constructor with explicit shift size and name override */ 
CEVPowerPointwise::CEVPowerPointwise(const string& overrideName,
                             double     shiftSize):
    VectorShift(TYPE, overrideName, shiftSize){
}

/** for reflection */
CEVPowerPointwise::CEVPowerPointwise(): 
    VectorShift(TYPE, NAME){
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double CEVPowerPointwise::divisor() const{
    // add check for 0 shift size - somewhere 
    return 100.0 * getShiftSize();
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP CEVPowerPointwise::shiftInterface() const{
    return IShift::TYPE;
}
 
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool CEVPowerPointwise::nameMatches(const OutputName&         name,
                        IObjectConstSP          obj){
    // cast obj to CEVPowerPointwise::Shift and then invoke name method
    const IShift& CEVPowerObj = 
        dynamic_cast<const IShift&>(*obj);
    return name.equals(CEVPowerObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's Shift interface */
void CEVPowerPointwise::appendName(OutputNameArray&          namesList,
                               IObjectConstSP          obj){
    // cast obj to CEVPowerPointwise::Shift and then invoke name method
    const IShift& CEVPowerObj = 
        dynamic_cast<const IShift&>(*obj);
    OutputNameSP outputName(new OutputName(CEVPowerObj.sensName(this)));
    namesList.push_back(outputName);
}

bool CEVPowerPointwise::shift(IObjectSP obj) {
    // cast obj to CEVPowerPointwise::Shift and then invoke shift method
    IShift& CEVPowerObj = 
        dynamic_cast<IShift&>(*obj);
    return CEVPowerObj.sensShift(this);
}

/** The supplied object is queried for the expiries array needed
    for doing CEVPower pointwise and this array is returned. The supplied
    object must implement the CEVPowerPointwise.Shift interface */
IObjectConstSP CEVPowerPointwise::qualifier(IObjectConstSP obj){
    // cast obj to CEVPowerPointwise::Shift and then invoke qualifier method
    const IShift& CEVPowerObj = 
        dynamic_cast<const IShift&>(*obj);
    return CEVPowerObj.sensExpiries(this);
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP CEVPowerPointwise::restorableShiftInterface() const{
    return IRestorableShift::TYPE;
}

void CEVPowerPointwise::restore(IObjectSP obj) {
    // cast obj to CEVPowerPointwise::Shift and then invoke restore method
    IRestorableShift& CEVPowerObj = 
        dynamic_cast<IRestorableShift&>(*obj);
    CEVPowerObj.sensRestore(this);
}

class CEVPowerPointwiseHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new CEVPowerPointwise(CEVPowerPointwise::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new CEVPowerPointwise(shiftSize);
        }
    };
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CEVPowerPointwise, clazz);
        SUPERCLASS(VectorShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultCEVPowerPointwise);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(CEVPowerPointwise::NAME, 
                                    new Factory(), 
                                    new CEVPowerPointwise(CEVPowerPointwise::DEFAULT_SHIFT),
                                    CEVPowerPointwise::IShift::TYPE);
    }

    static IObject* defaultCEVPowerPointwise(){
        return new CEVPowerPointwise();
    }
};

CClassConstSP const CEVPowerPointwise::TYPE = CClass::registerClassLoadMethod(
    "CEVPowerPointwise", typeid(CEVPowerPointwise), CEVPowerPointwiseHelper::load);

CClassConstSP const CEVPowerPointwise::IShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "CEVPowerPointwise::IShift", typeid(CEVPowerPointwise::IShift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(CEVPowerPointwise::IRestorableShift, clazz);
    EXTENDS(CEVPowerPointwise::IShift);
}
    
CClassConstSP const CEVPowerPointwise::IRestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "CEVPowerPointwise::IRestorableShift",
    typeid(CEVPowerPointwise::IRestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE
