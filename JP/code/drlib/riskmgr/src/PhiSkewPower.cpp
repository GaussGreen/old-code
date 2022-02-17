//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PhiSkewPower.cpp
//
//   Description : PHI_SKEWPOWER sensitivity
//
//   Author      : Oliver Brockhaus
//
//   Date        : 10 Oct 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/PhiSkewPower.hpp"
#include "edginc/Maths.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE

PhiSkewPower::Shift::~Shift(){} // empty
PhiSkewPower::RestorableShift::~RestorableShift(){} // empty

/** Sens Control for PHI_SKEW       */
const string PhiSkewPower::NAME = "PHI_SKEWPOWER";
const double PhiSkewPower::DEFAULT_SHIFT = 0.01;

/** constructor with explicit shift size */
PhiSkewPower::PhiSkewPower(double shiftSize):
    ScalarShift(TYPE, NAME, shiftSize){
}

/** for reflection */
PhiSkewPower::PhiSkewPower(): 
    ScalarShift(TYPE, NAME){
    validatePop2Object();
}

/** Validation */
void PhiSkewPower::validatePop2Object(){
    static const string method("PhiSkewPower::validatePop2Object");
    double shiftSize = getShiftSize();
    if ( (shiftSize > 0.5) || (shiftSize < 0.0) ){
        throw ModelException( method, "PhiSkewPower "+Format::toString(shiftSize)+
                              "must be between 0.0 and 0.5.");
    }
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double PhiSkewPower::divisor() const{
    static const char routine[] = "PhiSkewPower::divisor";
    double initVal;
    double shiftSize;
    double divisor;
    try{
        // find initial value of correlation skew - stored in sens control
        try{
            initVal = getInitialValue();
        } catch (ModelException& e){
            e.addMsg("Initial value of correlation skew has not been set"
                     " in the PHI_SKEWPOWER calculation");
            throw e;
        }
        shiftSize = getShiftSize();
        if (Maths::isZero(shiftSize))
        {
            throw ModelException(routine, "Shift size is zero");
        }
        // report a 0.1 move in correlation skew power; shift towards 0.5
        divisor = 10.0 * shiftSize * (initVal>0.5 ? -1.0 : 1.0);

    } catch (ModelException& e){
        e.addMsg(routine);
        throw e;
    }
    return divisor;
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP PhiSkewPower::shiftInterface() const{
    return Shift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP PhiSkewPower::restorableShiftInterface() const{
    return RestorableShift::TYPE;
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity.  The object must implement the
    VegaParallel.Shift interface */
bool PhiSkewPower::nameMatches(const OutputName&         name,
                        IObjectConstSP          obj){
    // cast obj to PhiSkewPower::Shift and then invoke name method
    const Shift& PhiSkewPowerObj =
        dynamic_cast<const Shift&>(*obj);
    return name.equals(PhiSkewPowerObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to
    this sensitivity to the supplied list */
void PhiSkewPower::appendName(OutputNameArray&          namesList,
                       IObjectConstSP          obj){
    // cast obj to PhiSkewPower::Shift and then invoke name method
    const Shift& PhiSkewPowerObj =
        dynamic_cast<const Shift&>(*obj);
    OutputNameSP outputName(new OutputName(PhiSkewPowerObj.sensName(this)));
    namesList.push_back(outputName);
}

bool PhiSkewPower::shift(IObjectSP obj) {
    // cast obj to PhiSkewPower::Shift and then invoke shift method
    Shift& PhiSkewPowerObj =
        dynamic_cast<Shift&>(*obj);
    return PhiSkewPowerObj.sensShift(this);
}

void PhiSkewPower::restore(IObjectSP obj) {
    // cast obj to PhiSkewPower::Shift and then invoke restore method
    RestorableShift& PhiSkewPowerObj =
        dynamic_cast<RestorableShift&>(*obj);
    PhiSkewPowerObj.sensRestore(this);
}

class PhiSkewPowerHelper{
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new PhiSkewPower(PhiSkewPower::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new PhiSkewPower(shiftSize);
        }
    };
    
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(PhiSkewPower, clazz);
        SUPERCLASS(ScalarShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultPhiSkewPower);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(PhiSkewPower::NAME, 
                                    new Factory(), 
                                    new PhiSkewPower(PhiSkewPower::DEFAULT_SHIFT),
                                    PhiSkewPower::Shift::TYPE);
    }

    static IObject* defaultPhiSkewPower(){
        return new PhiSkewPower();
    }
};

CClassConstSP const PhiSkewPower::TYPE = CClass::registerClassLoadMethod(
    "PhiSkewPower", typeid(PhiSkewPower), PhiSkewPowerHelper::load);

CClassConstSP const PhiSkewPower::Shift::TYPE =
CClass::registerInterfaceLoadMethod(
    "PhiSkewPower::Shift", typeid(PhiSkewPower::Shift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(PhiSkewPower::RestorableShift, clazz);
    EXTENDS(PhiSkewPower::Shift);
}
    
CClassConstSP const PhiSkewPower::RestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "PhiSkewPower::RestorableShift", typeid(PhiSkewPower::RestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE






