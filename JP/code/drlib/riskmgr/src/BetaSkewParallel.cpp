//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids Derivatives Research
//
//   Filename    : BetaSkewParallel.hpp
//
//   Description : Pointwise Shift on a beta skew matrix
//
//   Author      : Gordon Stephens
//
//   Date        : 23 March 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/BetaSkewParallel.hpp"
#include "edginc/Maths.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE
BetaSkewParallel::IShift::~IShift(){} // empty
BetaSkewParallel::IRestorableShift::~IRestorableShift(){} // empty

const double ONE_BASIS_POINT = 0.0001;

/** Sens Control for rho parallel */
const string BetaSkewParallel::NAME = "BETA_SKEW_PARALLEL";
const double BetaSkewParallel::DEFAULT_SHIFT = ONE_BASIS_POINT;

/** constructor with explicit shift size */
BetaSkewParallel::BetaSkewParallel(double shiftSize):
    ScalarShift(TYPE, NAME, shiftSize){
}

/** constructor with explicit shift size */
BetaSkewParallel::BetaSkewParallel(double     shiftSize,
             IModel*    model,
             Control*   control):
    ScalarShift(TYPE, NAME, shiftSize) {
    this->algorithm = model;
    this->control   = control;
}

/** for reflection */
BetaSkewParallel::BetaSkewParallel(): 
    ScalarShift(TYPE, NAME){
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double BetaSkewParallel::divisor() const{
    static const string method = "BetaSkewParallel::divisor";
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
CClassConstSP BetaSkewParallel::shiftInterface() const{
    return IShift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP BetaSkewParallel::restorableShiftInterface() const{
    return IRestorableShift::TYPE;
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement the
    Shift interface */
bool BetaSkewParallel::nameMatches(const OutputName&         name,
                               IObjectConstSP            obj){
    // cast obj to BetaSkewParallel::IShift and then invoke name method
    const IShift& skewParallelObj =
        dynamic_cast<const IShift&>(*obj);
    return name.equals(skewParallelObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement the
    Shift interface */
void BetaSkewParallel::appendName(OutputNameArray&          namesList,
                              IObjectConstSP            obj){
    // cast obj to BetaSkewParallel::IShift and then invoke name method
    const IShift& skewParallelObj =
        dynamic_cast<const IShift&>(*obj);
    OutputNameSP outputName(new OutputName(skewParallelObj.sensName(this)));
    namesList.push_back(outputName);
}

bool BetaSkewParallel::shift(IObjectSP obj) {
    // cast obj to RhoParallel::Shift and then invoke shift method
    IShift& skewParallelObj =
        dynamic_cast<IShift&>(*obj);
    return skewParallelObj.sensShift(this);
}

    
void BetaSkewParallel::restore(IObjectSP obj) {
    // cast obj to BetaSkewParallel::RestorableShift and then invoke restore method
    IRestorableShift& skewParallelObj =
        dynamic_cast<IRestorableShift&>(*obj);
    skewParallelObj.sensRestore(this);
}

class BetaSkewParallelHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new BetaSkewParallel(BetaSkewParallel::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new BetaSkewParallel(shiftSize);
        }
    };
        
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(BetaSkewParallel, clazz);
        SUPERCLASS(ScalarShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultBetaSkewParallel);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(BetaSkewParallel::NAME, 
                                    new Factory(), 
                                    new BetaSkewParallel(BetaSkewParallel::DEFAULT_SHIFT),
                                    BetaSkewParallel::IShift::TYPE);
    }

    static IObject* defaultBetaSkewParallel(){
        return new BetaSkewParallel();
    }
};

CClassConstSP const BetaSkewParallel::TYPE = CClass::registerClassLoadMethod(
    "BetaSkewParallel", typeid(BetaSkewParallel), BetaSkewParallelHelper::load);

CClassConstSP const BetaSkewParallel::IShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "BetaSkewParallel::IShift", typeid(BetaSkewParallel::IShift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(BetaSkewParallel::IRestorableShift, clazz);
    EXTENDS(BetaSkewParallel::IShift);
}
    
CClassConstSP const BetaSkewParallel::IRestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "BetaSkewParallel::IRestorableShift", typeid(BetaSkewParallel::IRestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE
