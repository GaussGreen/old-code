
#include "edginc/config.hpp"
#include "edginc/CreditSpreadVegaParallel.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE
CreditSpreadVegaParallel::Shift::~Shift(){} // empty
CreditSpreadVegaParallel::RestorableShift::~RestorableShift(){} // empty

/** Sens Control for vega parallel */
const string CreditSpreadVegaParallel::NAME = "CREDIT_SPREAD_VEGA_PARALLEL";
const double CreditSpreadVegaParallel::DEFAULT_SHIFT = 0.001;

/** constructor with explicit shift size */
CreditSpreadVegaParallel::CreditSpreadVegaParallel(double     shiftSize):
    ScalarShift(TYPE, NAME, shiftSize){
}

/** for reflection */
CreditSpreadVegaParallel::CreditSpreadVegaParallel(): 
    ScalarShift(TYPE, NAME){
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double CreditSpreadVegaParallel::divisor() const{
    // add check for 0 shift size - somewhere 
    return 100.0 * getShiftSize();
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP CreditSpreadVegaParallel::shiftInterface() const{
    return Shift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP CreditSpreadVegaParallel::restorableShiftInterface() const{
    return RestorableShift::TYPE;
}

/** Returns true if the supplied object matches the supplied name
    for this sensitivity.  The object must implement the
    CreditSpreadVegaParallel.Shift interface */
bool CreditSpreadVegaParallel::nameMatches(const OutputName&         name,
                                    IObjectConstSP            obj){
    // cast obj to CreditSpreadVegaParallel::Shift and then invoke name method
    const Shift& vegaObj =
        dynamic_cast<const Shift&>(*obj);
    return name.equals(vegaObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to
    this sensitivity to the supplied list */
void CreditSpreadVegaParallel::appendName(OutputNameArray&          namesList,
                                   IObjectConstSP            obj){
    // cast obj to CreditSpreadVegaParallel::Shift and then invoke name method
    const Shift& vegaObj =
        dynamic_cast<const Shift&>(*obj);
    OutputNameSP outputName(new OutputName(vegaObj.sensName(this)));
    namesList.push_back(outputName);
}

bool CreditSpreadVegaParallel::shift(IObjectSP obj) {
    // cast obj to CreditSpreadVegaParallel::Shift and then invoke shift method
    Shift& vegaObj =
        dynamic_cast<Shift&>(*obj);
    return vegaObj.sensShift(this);
}

    
void CreditSpreadVegaParallel::restore(IObjectSP obj) {
    // cast obj to CreditSpreadVegaParallel::Shift and then invoke restore method
    RestorableShift& vegaObj =
        dynamic_cast<RestorableShift&>(*obj);
    vegaObj.sensRestore(this);
}

class CreditSpreadVegaParallelHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new CreditSpreadVegaParallel(CreditSpreadVegaParallel::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new CreditSpreadVegaParallel(shiftSize);
        }
    };
        
public:
    /** Invoked when this class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CreditSpreadVegaParallel, clazz);
        SUPERCLASS(ScalarShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultCreditSpreadVegaParallel);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(CreditSpreadVegaParallel::NAME, 
                                    new Factory(), 
                                    new CreditSpreadVegaParallel(CreditSpreadVegaParallel::DEFAULT_SHIFT),
                                    CreditSpreadVegaParallel::Shift::TYPE);
    }

    static IObject* defaultCreditSpreadVegaParallel(){
        return new CreditSpreadVegaParallel();
    }
};

CClassConstSP const CreditSpreadVegaParallel::TYPE = CClass::registerClassLoadMethod(
    "CreditSpreadVegaParallel", typeid(CreditSpreadVegaParallel), CreditSpreadVegaParallelHelper::load);

CClassConstSP const CreditSpreadVegaParallel::Shift::TYPE =
CClass::registerInterfaceLoadMethod(
    "CreditSpreadVegaParallel::Shift", typeid(CreditSpreadVegaParallel::Shift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(CreditSpreadVegaParallel::RestorableShift, clazz);
    EXTENDS(CreditSpreadVegaParallel::Shift);
}
    
CClassConstSP const CreditSpreadVegaParallel::RestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "CreditSpreadVegaParallel::RestorableShift", typeid(CreditSpreadVegaParallel::RestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE
