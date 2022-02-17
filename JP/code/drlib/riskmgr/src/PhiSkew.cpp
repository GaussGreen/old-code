//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PhiSkew.cpp
//
//   Description : PHI_SKEW sensitivity
//
//   Author      : Oliver Brockhaus
//
//   Date        : 10 Oct 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/PhiSkew.hpp"
#include "edginc/Maths.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE

PhiSkew::Shift::~Shift(){} // empty
PhiSkew::RestorableShift::~RestorableShift(){} // empty

/** Sens Control for PHI_SKEW       */
const string PhiSkew::NAME = "PHI_SKEW";
const double PhiSkew::DEFAULT_SHIFT = 0.01;

/** constructor with explicit shift size */
PhiSkew::PhiSkew(double shiftSize):
    ScalarShift(TYPE, NAME, shiftSize){
    validatePop2Object();
}

/** for reflection */
PhiSkew::PhiSkew(): 
    ScalarShift(TYPE, NAME){
}

/** Validation */
void PhiSkew::validatePop2Object(){
    static const string method("PhiSkew::validatePop2Object");
    double shiftSize = getShiftSize();
    if ( (shiftSize > 0.5) || (shiftSize < 0.0) ){
        throw ModelException( method, "PhiSkew "+Format::toString(shiftSize)+
                              "must be between 0.0 and 0.5.");
    }
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double PhiSkew::divisor() const{
    static const char routine[] = "PhiSkew::divisor";
    double initVal;
    double shiftSize;
    double divisor;
    try{
        // find initial value of correlation skew - stored in sens control
        try{
            initVal = getInitialValue();
        } catch (ModelException& e){
            e.addMsg("Initial value of correlation skew has not been set"
                     " in the PHI_SKEW calculation");
            throw e;
        }
        // then just scale by shift size
        shiftSize = getShiftSize();
        if (Maths::isZero(shiftSize))
        {
            throw ModelException(routine, "Shift size is zero");
        }
        // report a 0.1 move in correlation skew (in [0,1]); shift towards 0.5
        divisor = 10.0 * shiftSize * (initVal>0.5 ? -1.0 : 1.0);

    } catch (ModelException& e){
        e.addMsg(routine);
        throw e;
    }
    return divisor;
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP PhiSkew::shiftInterface() const{
    return Shift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP PhiSkew::restorableShiftInterface() const{
    return RestorableShift::TYPE;
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity.  The object must implement the
    VegaParallel.Shift interface */
bool PhiSkew::nameMatches(const OutputName&         name,
                        IObjectConstSP          obj){
    // cast obj to PhiSkew::Shift and then invoke name method
    const Shift& PhiSkewObj =
        dynamic_cast<const Shift&>(*obj);
    return name.equals(PhiSkewObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to
    this sensitivity to the supplied list */
void PhiSkew::appendName(OutputNameArray&          namesList,
                       IObjectConstSP          obj){
    // cast obj to PhiSkew::Shift and then invoke name method
    const Shift& PhiSkewObj =
        dynamic_cast<const Shift&>(*obj);
    OutputNameSP outputName(new OutputName(PhiSkewObj.sensName(this)));
    namesList.push_back(outputName);
}

bool PhiSkew::shift(IObjectSP obj) {
    // cast obj to PhiSkew::Shift and then invoke shift method
    Shift& PhiSkewObj =
        dynamic_cast<Shift&>(*obj);
    return PhiSkewObj.sensShift(this);
}

void PhiSkew::restore(IObjectSP obj) {
    // cast obj to PhiSkew::Shift and then invoke restore method
    RestorableShift& PhiSkewObj =
        dynamic_cast<RestorableShift&>(*obj);
    PhiSkewObj.sensRestore(this);
}

class PhiSkewHelper{
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new PhiSkew(PhiSkew::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new PhiSkew(shiftSize);
        }
    };
    
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(PhiSkew, clazz);
        SUPERCLASS(ScalarShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultPhiSkew);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(PhiSkew::NAME, 
                                    new Factory(), 
                                    new PhiSkew(PhiSkew::DEFAULT_SHIFT),
                                    PhiSkew::Shift::TYPE);
    }

    static IObject* defaultPhiSkew(){
        return new PhiSkew();
    }
};

CClassConstSP const PhiSkew::TYPE = CClass::registerClassLoadMethod(
    "PhiSkew", typeid(PhiSkew), PhiSkewHelper::load);

CClassConstSP const PhiSkew::Shift::TYPE =
CClass::registerInterfaceLoadMethod(
    "PhiSkew::Shift", typeid(PhiSkew::Shift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(PhiSkew::RestorableShift, clazz);
    EXTENDS(PhiSkew::Shift);
}
    
CClassConstSP const PhiSkew::RestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "PhiSkew::RestorableShift", typeid(PhiSkew::RestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE






