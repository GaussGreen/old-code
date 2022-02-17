//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FXVega.cpp
//
//   Description : Controls calculation of FX_VEGA
//
//   Author      : Mark A Robson
//
//   Date        : 15 Mar 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/FXVega.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE
FXVega::IShift::~IShift(){} // empty
FXVega::IRestorableShift::~IRestorableShift(){} // empty

/** Sens Control for vega parallel */
const string FXVega::NAME = "FX_VEGA";
const double FXVega::DEFAULT_SHIFT = 0.001;

/** constructor with explicit shift size */
FXVega::FXVega(double     shiftSize):
    ScalarShift(TYPE, NAME, shiftSize){
}

FXVega::FXVega(CClassConstSP clazz, const string& sensName): 
    ScalarShift(clazz, sensName) {}

/** for reflection */
FXVega::FXVega(): 
    ScalarShift(TYPE, NAME){
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double FXVega::divisor() const{
    // add check for 0 shift size - somewhere 
    return 100.0 * getShiftSize();
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP FXVega::shiftInterface() const{
    return IShift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP FXVega::restorableShiftInterface() const{
    return IRestorableShift::TYPE;
}

/** Returns true if the supplied object matches the supplied name
    for this sensitivity.  The object must implement the
    FXVega.Shift interface */
bool FXVega::nameMatches(const OutputName&         name,
                               IObjectConstSP          obj){
    // cast obj to FXVega::Shift and then invoke name method
    const IShift& vegaObj =
        dynamic_cast<const IShift&>(*obj);
    return name.equals(vegaObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to
    this sensitivity to the supplied list */
void FXVega::appendName(OutputNameArray&          namesList,
                        IObjectConstSP          obj){
    // cast obj to FXVega::Shift and then invoke name method
    const IShift& vegaObj =
        dynamic_cast<const IShift&>(*obj);
    OutputNameSP outputName(new OutputName(vegaObj.sensName(this)));
    namesList.push_back(outputName);
}

bool FXVega::shift(IObjectSP obj) {
    // cast obj to FXVega::Shift and then invoke shift method
    IShift& vegaObj =
        dynamic_cast<IShift&>(*obj);
    return vegaObj.sensShift(this);
}

void FXVega::restore(IObjectSP obj) {
    // cast obj to FXVega::Shift and then invoke restore method
    IRestorableShift& vegaObj =
        dynamic_cast<IRestorableShift&>(*obj);
    vegaObj.sensRestore(this);
}

// so we get some naming consistency with everything else
class FXVegaParallel: public FXVega {
public:
    static CClassConstSP const TYPE;
    const static string NAME;
    
private:
    FXVegaParallel():FXVega(TYPE, NAME) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(FXVegaParallel, clazz);
        SUPERCLASS(FXVega);
        EMPTY_SHELL_METHOD(defaultFXVegaParallel);
        clazz->setPublic(); // make visible to EAS/spreadsheet
        clazz->setDescription("simple parallel shift - same as FXVega");
    }

    static IObject* defaultFXVegaParallel(){
        return new FXVegaParallel();
    }
};

const string FXVegaParallel::NAME = "FX_VEGA_PARALLEL";

CClassConstSP const FXVegaParallel::TYPE = CClass::registerClassLoadMethod(
    "FXVegaParallel", typeid(FXVegaParallel), load);


class FXVegaHelper{
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new FXVega(FXVega::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new FXVega(shiftSize);
        }
    };
        
public:
    /** Invoked when SmallClass is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(FXVega, clazz);
        SUPERCLASS(ScalarShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultFXVega);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(FXVega::NAME, 
                                    new Factory(), 
                                    new FXVega(FXVega::DEFAULT_SHIFT),
                                    FXVega::IShift::TYPE);
        SensitivityFactory::addSens(FXVegaParallel::NAME, 
                                    new Factory(), 
                                    new FXVega(FXVega::DEFAULT_SHIFT),
                                    FXVegaParallel::IShift::TYPE);
    }

    static IObject* defaultFXVega(){
        return new FXVega();
    }
};

CClassConstSP const FXVega::TYPE = CClass::registerClassLoadMethod(
    "FXVega", typeid(FXVega), FXVegaHelper::load);

CClassConstSP const FXVega::IShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "FXVega::IShift", typeid(FXVega::IShift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(FXVega::IRestorableShift, clazz);
    EXTENDS(FXVega::IShift);
}
    
CClassConstSP const FXVega::IRestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "FXVega::IRestorableShift", typeid(FXVega::IRestorableShift), 
    restorableShiftLoad);


DRLIB_END_NAMESPACE
