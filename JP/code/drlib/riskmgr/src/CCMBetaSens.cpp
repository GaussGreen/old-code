//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CCMBetaSens.cpp
//
//   Description : Sensitivity corresponding to shifting CCM Beta parameter
//
//   Author      : Antoine Gregoire
//
//   Date        : October 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CCMBetaSens.hpp"
#include "edginc/Maths.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE

CCMBetaSens::IShift::~IShift(){} // empty
CCMBetaSens::IRestorableShift::~IRestorableShift(){} // empty

/** Sens Control for beta */
const string CCMBetaSens::NAME = "CCM_BETA_SENS";
const double CCMBetaSens::DEFAULT_SHIFT = 0.001;//TODO check
 
CCMBetaSens::~CCMBetaSens(){}

/** constructor with explicit shift size */
CCMBetaSens::CCMBetaSens(double     shiftSize):
    ScalarShift(TYPE, NAME, shiftSize){}

/** for reflection */
CCMBetaSens::CCMBetaSens(): ScalarShift(TYPE, NAME){}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */

double CCMBetaSens::divisor() const {
    double old;
    try {
        old = getInitialValue();
    }
    catch (exception& e) {
        throw ModelException(e, "CCMBetaSens::divisor");
    }
    double relShift = getShiftSize();
    if (Maths::isZero(relShift)) {
        throw ModelException("CCMBetaSens::divisor", "Shift size is zero");
    }
    double it = relShift * (1 - old);
    return Maths::isZero(it) ? 1. : it;
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP CCMBetaSens::shiftInterface() const{
    return IShift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP CCMBetaSens::restorableShiftInterface() const{
    return IRestorableShift::TYPE;
}

/** Returns true if the supplied object matches the supplied name
    for this sensitivity.  The object must implement the
    CCMBetaSens.Shift interface */
bool CCMBetaSens::nameMatches(const OutputName&         name,
                              IObjectConstSP            obj){
    // cast obj to CCMBetaSens::Shift and then invoke name method
    const IShift& betaObj = dynamic_cast<const IShift&>(*obj);
    return name.equals(betaObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to
    this sensitivity to the supplied list */
void CCMBetaSens::appendName(OutputNameArray&          namesList,
                             IObjectConstSP            obj){
    // cast obj to CCMBetaSens::Shift and then invoke name method

  
    const IShift& betaObj = dynamic_cast<const IShift&>(*obj);
    OutputNameSP outputName(new OutputName(betaObj.sensName(this)));
    namesList.push_back(outputName);
}

bool CCMBetaSens::shift(IObjectSP obj) {
    // cast obj to CCMBetaSens::Shift and then invoke shift method
    IShift& betaObj = dynamic_cast<IShift&>(*obj);
    return betaObj.sensShift(this);
}

    
void CCMBetaSens::restore(IObjectSP obj) {
    // cast obj to CCMBetaSens::Shift and then invoke restore method
    IRestorableShift& betaObj = dynamic_cast<IRestorableShift&>(*obj);
    betaObj.sensRestore(this);
}

class CCMBetaSensHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public virtual SensitivityFactory::IDefault,
                   public virtual SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new CCMBetaSens(CCMBetaSens::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new CCMBetaSens(shiftSize);
        }
    };
        
public:
    /** Invoked when this class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CCMBetaSens, clazz);
        SUPERCLASS(ScalarShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultCCMBetaSens);
        // no fields (shift size is in parent)
        // register how to build our sensitivity
        SensitivityFactory::addSens(CCMBetaSens::NAME, 
                                    new Factory(), 
                                    new CCMBetaSens(CCMBetaSens::DEFAULT_SHIFT),
                                    CCMBetaSens::IShift::TYPE);
    }

    static IObject* defaultCCMBetaSens(){
        return new CCMBetaSens();
    }
};

CClassConstSP const CCMBetaSens::TYPE = CClass::registerClassLoadMethod(
    "CCMBetaSens", typeid(CCMBetaSens), CCMBetaSensHelper::load);

CClassConstSP const CCMBetaSens::IShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "CCMBetaSens::Shift", typeid(CCMBetaSens::IShift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(CCMBetaSens::IRestorableShift, clazz);
    EXTENDS(CCMBetaSens::IShift);
}
    
CClassConstSP const CCMBetaSens::IRestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "CCMBetaSens::RestorableShift", typeid(CCMBetaSens::IRestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE
