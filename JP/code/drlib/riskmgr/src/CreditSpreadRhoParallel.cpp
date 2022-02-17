//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CreditSpreadRhoParallel.cpp
//
//   Description : spread rho parallel sensitivity
//
//   Author      : André Segger
//
//   Date        : 16 April 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CreditSpreadRhoParallel.hpp"
#include "edginc/Maths.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/Results.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/TweakGroup.hpp"

DRLIB_BEGIN_NAMESPACE
typedef smartPtr<CreditSpreadRhoParallel::RestorableShift> CreditSpreadRhoParallelRestorableShiftSP;

// not sure why we have to have these - MSVC gets upset if they're not in the
// header
CreditSpreadRhoParallel::Shift::~Shift(){} // empty
CreditSpreadRhoParallel::RestorableShift::~RestorableShift(){} // empty

const double ONE_BASIS_POINT = 0.0001;

/** Sens Control for rho parallel */
const string CreditSpreadRhoParallel::NAME = "CREDIT_SPREAD_RHO_PARALLEL";
const double CreditSpreadRhoParallel::DEFAULT_SHIFT = ONE_BASIS_POINT;

/** constructor with explicit shift size */
CreditSpreadRhoParallel::CreditSpreadRhoParallel(double shiftSize):
    ScalarShift(TYPE, NAME, shiftSize){
}

/** for reflection */
CreditSpreadRhoParallel::CreditSpreadRhoParallel(): 
    ScalarShift(TYPE, NAME){
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double CreditSpreadRhoParallel::divisor() const{
    static const string method = "CreditSpreadRhoParallel::divisor";
    double shiftSize;
    try{
        // just scale the shift size for a 1bp move
        shiftSize = getShiftSize();
        if (Maths::isZero(shiftSize)){
            throw ModelException(method, "Shift size is zero");
        }
        return (shiftSize/ONE_BASIS_POINT);
    } 
    catch (ModelException& e) {
        throw ModelException(&e, method);
    }
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP CreditSpreadRhoParallel::shiftInterface() const{
    return Shift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP CreditSpreadRhoParallel::restorableShiftInterface() const{
    return RestorableShift::TYPE;
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement the
    Shift interface */
bool CreditSpreadRhoParallel::nameMatches(const OutputName&         name,
                                          IObjectConstSP          obj) {
    // cast obj to RhoParallel::Shift and then invoke name method
    const Shift& rhoParallelObj = 
        dynamic_cast<const Shift&>(*obj);
    return name.equals(rhoParallelObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement the
    Shift interface */
void CreditSpreadRhoParallel::appendName(OutputNameArray&          namesList,
                                         IObjectConstSP          obj){
    // cast obj to RhoParallel::Shift and then invoke name method
    const Shift& rhoParallelObj = 
        dynamic_cast<const Shift&>(*obj);
    OutputNameSP outputName(new OutputName(rhoParallelObj.sensName(this)));
    namesList.push_back(outputName);
}

bool CreditSpreadRhoParallel::shift(IObjectSP obj) {
    // cast obj to RhoParallel::Shift and then invoke shift method
    Shift& rhoParallelObj = 
        dynamic_cast<Shift&>(*obj);
    return rhoParallelObj.sensShift(this);
}

void CreditSpreadRhoParallel::restore(IObjectSP obj) {
    // cast obj to RhoParallel::Shift and then invoke restore method
    RestorableShift& rhoParallelObj = 
        dynamic_cast<RestorableShift&>(*obj);
    rhoParallelObj.sensRestore(this);
}

/** implements a one sided scalar derivative for each instance of the
    market data which is sensitive to this SensControl */
void CreditSpreadRhoParallel::calculate(TweakGroup*  tweakGroup,
                                        CResults*    results){
    try {
        // get list of names to calculate result for
        OutputNameArrayConstSP names(this->names(tweakGroup, results));

        if (names->empty()) {
            results->storeNotApplicable(this);            
        } else {
            double origPrice;

            for (int idx = 0; idx < names->size(); idx++){
                // store what we want to shift
                setMarketDataName((*names)[idx]);

                // We need to compute a base price for pricing off credit spreads before tweaking
                // We achieve this by shifting and pricing with 0 shiftSize
                double restoreShiftSize = getShiftSize();
                setShiftSize(0.0);

                try {
                    double dummy = 0.0;     // Not used
                    origPrice = shiftAndPrice(tweakGroup, dummy, true);
                } 
                catch (exception& e){
                    setShiftSize(restoreShiftSize);
                    throw ModelException(e);
                }
                setShiftSize(restoreShiftSize);

                /* skip over where result has been calculated already */
                if (!results->exists(this)){
                    try {
                        // calculate sens
                        double shiftedPrice = shiftAndPrice(tweakGroup, origPrice);
                        double divisor = this->divisor();
                        double firstDeriv = (shiftedPrice - origPrice)/divisor;

                        // and store it
                        results->storeScalarGreek(firstDeriv, this);
                    }
                    catch (exception& e) {
                        results->storeGreek(IObjectSP(new Untweakable(e)), this);
                    }
                }
            }
        }
    } catch (exception& e){
        results->storeGreek(IObjectSP(new Untweakable(e)), getSensOutputName(),
                            OutputNameSP(new OutputName("")));
    }
}

class CreditSpreadRhoParallelHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new CreditSpreadRhoParallel(CreditSpreadRhoParallel::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new CreditSpreadRhoParallel(shiftSize);
        }
    };
        
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CreditSpreadRhoParallel, clazz);
        SUPERCLASS(ScalarShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultRhoParallel);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(CreditSpreadRhoParallel::NAME, 
                                    new Factory(), 
                                    new CreditSpreadRhoParallel(CreditSpreadRhoParallel::DEFAULT_SHIFT),
                                    CreditSpreadRhoParallel::Shift::TYPE);
    }

    static IObject* defaultRhoParallel(){
        return new CreditSpreadRhoParallel();
    }
};

CClassConstSP const CreditSpreadRhoParallel::TYPE = CClass::registerClassLoadMethod(
    "CreditSpreadRhoParallel", typeid(CreditSpreadRhoParallel), CreditSpreadRhoParallelHelper::load);

CClassConstSP const CreditSpreadRhoParallel::Shift::TYPE =
CClass::registerInterfaceLoadMethod(
    "CreditSpreadRhoParallel::Shift", typeid(CreditSpreadRhoParallel::Shift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(CreditSpreadRhoParallel::RestorableShift, clazz);
    EXTENDS(CreditSpreadRhoParallel::Shift);
}
    
CClassConstSP const CreditSpreadRhoParallel::RestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "CreditSpreadRhoParallel::RestorableShift", typeid(CreditSpreadRhoParallel::RestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE
