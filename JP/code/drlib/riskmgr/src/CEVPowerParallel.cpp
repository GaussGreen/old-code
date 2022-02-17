//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CEVPowerParallel.cpp
//
//   Description : for CEVPower parallel tweak
//
//   Date        : 25 Feb 2002
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/CEVPowerParallel.hpp"

DRLIB_BEGIN_NAMESPACE

// not sure why we have to have these - MSVC gets upset if they're not in the
// header
CEVPowerParallel::Shift::~Shift(){} // empty
CEVPowerParallel::RestorableShift::~RestorableShift(){} // empty

/** Sens Control for CEVPower parallel */
const string CEVPowerParallel::NAME = "CEVPOWER_PARALLEL";
const double CEVPowerParallel::DEFAULT_SHIFT = 0.1;

/** constructor with explicit shift size */
CEVPowerParallel::CEVPowerParallel(double     shiftSize):
    ScalarShift(TYPE, NAME, shiftSize){
}

/** for reflection */
CEVPowerParallel::CEVPowerParallel(): 
    ScalarShift(TYPE, NAME){
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double CEVPowerParallel::divisor() const{
    // add check for 0 shift size - somewhere 
    return 100.0 * getShiftSize();
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP CEVPowerParallel::shiftInterface() const{
    return Shift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP CEVPowerParallel::restorableShiftInterface() const{
    return RestorableShift::TYPE;
}

/** Returns true if the supplied object matches the supplied name
    for this sensitivity.  The object must implement the
    CEVPowerParallel.Shift interface */
bool CEVPowerParallel::nameMatches(const OutputName&         name,
                               IObjectConstSP          obj){
    // cast obj to CEVPowerParallel::Shift and then invoke name method
    const Shift& CEVPowerObj = 
        dynamic_cast<const Shift&>(*obj);
    return name.equals(CEVPowerObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to
    this sensitivity to the supplied list */
void CEVPowerParallel::appendName(OutputNameArray&          namesList,
                              IObjectConstSP          obj){
    // cast obj to CEVPowerParallel::Shift and then invoke name method
    const Shift& CEVPowerObj = 
        dynamic_cast<const Shift&>(*obj);
    OutputNameSP outputName(new OutputName(CEVPowerObj.sensName(this)));
    namesList.push_back(outputName);
}

bool CEVPowerParallel::shift(IObjectSP obj) {
    // cast obj to CEVPowerParallel::Shift and then invoke shift method
    Shift& CEVPowerObj = 
        dynamic_cast<Shift&>(*obj);
    return CEVPowerObj.sensShift(this);
}

void CEVPowerParallel::restore(IObjectSP obj) {
    // cast obj to CEVPowerParallel::Shift and then invoke restore method
    RestorableShift& CEVPowerObj = 
        dynamic_cast<RestorableShift&>(*obj);
    CEVPowerObj.sensRestore(this);
}

class CEVPowerParallelHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new CEVPowerParallel(CEVPowerParallel::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new CEVPowerParallel(shiftSize);
        }
    };
        
public:
    /** Invoked when this class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CEVPowerParallel, clazz);
        SUPERCLASS(ScalarShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultCEVPowerParallel);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(CEVPowerParallel::NAME, 
                                    new Factory(), 
                                    new CEVPowerParallel(CEVPowerParallel::DEFAULT_SHIFT),
                                    CEVPowerParallel::Shift::TYPE);
    }

    static IObject* defaultCEVPowerParallel(){
        return new CEVPowerParallel();
    }
};

CClassConstSP const CEVPowerParallel::TYPE = CClass::registerClassLoadMethod(
    "CEVPowerParallel", typeid(CEVPowerParallel), CEVPowerParallelHelper::load);

CClassConstSP const CEVPowerParallel::Shift::TYPE =
CClass::registerInterfaceLoadMethod(
    "CEVPowerParallel::Shift", typeid(CEVPowerParallel::Shift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(CEVPowerParallel::RestorableShift, clazz);
    EXTENDS(CEVPowerParallel::Shift);
}
    
CClassConstSP const CEVPowerParallel::RestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "CEVPowerParallel::RestorableShift", typeid(CEVPowerParallel::RestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE
