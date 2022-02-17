//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RootTimeVega.cpp
//
//   Description : Controls calculation of Root TIme Vega
//
//   Author      : Mark A Robson
//
//   Date        : 7 March 2001
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/RootTimeVega.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE
RootTimeVega::IShift::IShift(){} // empty
RootTimeVega::IShift::~IShift(){} // empty
RootTimeVega::IRestorableShift::IRestorableShift(){} // empty
RootTimeVega::IRestorableShift::~IRestorableShift(){} // empty

/** Sens Control for vega parallel */
const string RootTimeVega::NAME = "RT_VEGA";
const double RootTimeVega::DEFAULT_SHIFT = 0.001;

/** constructor with explicit shift size */
RootTimeVega::RootTimeVega(double     shiftSize):
    ScalarShift(TYPE, NAME, shiftSize){
}

/** for reflection */
RootTimeVega::RootTimeVega(): 
    ScalarShift(TYPE, NAME){
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double RootTimeVega::divisor() const{
    // add check for 0 shift size - somewhere 
    return 100.0 * getShiftSize();
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP RootTimeVega::shiftInterface() const{
    return IShift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP RootTimeVega::restorableShiftInterface() const{
    return IRestorableShift::TYPE;
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool RootTimeVega::nameMatches(const OutputName&         name,
                        IObjectConstSP          obj){
    // cast obj to RootTimeVega::Shift and then invoke name method
     const IShift& rootTimeVegaObj =
        dynamic_cast<const IShift&>(*obj);
     return name.equals(rootTimeVegaObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's Shift interface */
void RootTimeVega::appendName(OutputNameArray&          namesList,
                       IObjectConstSP          obj){
    // cast obj to RootTimeVega::Shift and then invoke name method
    const IShift& rootTimeVegaObj =
        dynamic_cast<const IShift&>(*obj);
    OutputNameSP outputName(new OutputName(rootTimeVegaObj.sensName(this)));
    namesList.push_back(outputName);
}

bool RootTimeVega::shift(IObjectSP obj) {
    // cast obj to RootTimeVega::IShift and then invoke shift method
    IShift& vegaObj =
        dynamic_cast<IShift&>(*obj);
    return vegaObj.sensShift(this);
}

bool RootTimeVega::restorableShift(IObjectConstSP& obj) const{
    return IRestorableShift::TYPE->isInstance(obj);
}
    
void RootTimeVega::restore(IObjectSP obj) {
    // cast obj to RootTimeVega::Shift and then invoke restore method
    IRestorableShift& vegaObj =
        dynamic_cast<IRestorableShift&>(*obj);
    vegaObj.sensRestore(this);
}

/** return shift/sqrt(time) */
double RootTimeVega::rtVegaShift(double years) const {
    if (years > 0.0)
    {
        return (getShiftSize()/sqrt(years));
    }
    else
    {   // for maturity = 0, don't shift
        return 0.0;
    }
}

class RootTimeVegaHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new RootTimeVega(RootTimeVega::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new RootTimeVega(shiftSize);
        }
    };
public:
    /** Invoked when SmallClass is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(RootTimeVega, clazz);
        SUPERCLASS(ScalarShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultRootTimeVega);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(RootTimeVega::NAME, 
                                    new Factory(), 
                                    new RootTimeVega(RootTimeVega::DEFAULT_SHIFT),
                                    RootTimeVega::IShift::TYPE);
    }

    static IObject* defaultRootTimeVega(){
        return new RootTimeVega();
    }
};

CClassConstSP const RootTimeVega::TYPE = CClass::registerClassLoadMethod(
    "RootTimeVega", typeid(RootTimeVega), RootTimeVegaHelper::load);

CClassConstSP const RootTimeVega::IShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "RootTimeVega::IShift", typeid(RootTimeVega::IShift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(RootTimeVega::IRestorableShift, clazz);
    EXTENDS(RootTimeVega::IShift);
}
    
CClassConstSP const RootTimeVega::IRestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "RootTimeVega::IRestorableShift", typeid(RootTimeVega::IRestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE
