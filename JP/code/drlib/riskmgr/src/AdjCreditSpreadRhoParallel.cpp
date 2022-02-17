//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : AdjCreditSpreadRhoParallel.cpp
//
//   Description : Adjusted credit spread rho parallel sensitivity
//
//   Author      : Milan Kovacevic
//
//   Date        : 20 November 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/AdjCreditSpreadRhoParallel.hpp"
#include "edginc/Maths.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/Results.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/Holiday.hpp"

DRLIB_BEGIN_NAMESPACE
typedef smartPtr<AdjCreditSpreadRhoParallel::RestorableShift> AdjCreditSpreadRhoParallelRestorableShiftSP;

// not sure why we have to have these - MSVC gets upset if they're not in the
// header
AdjCreditSpreadRhoParallel::Shift::~Shift(){} // empty
AdjCreditSpreadRhoParallel::RestorableShift::~RestorableShift(){} // empty

const double ONE_BASIS_POINT = 0.0001;

/** Sens Control for nbrho parallel */
const string AdjCreditSpreadRhoParallel::NAME = "ADJ_CREDIT_SPREAD_RHO_PARALLEL";
const double AdjCreditSpreadRhoParallel::DEFAULT_SHIFT = ONE_BASIS_POINT;

/** constructor with explicit shift size */
AdjCreditSpreadRhoParallel::AdjCreditSpreadRhoParallel(double shiftSize):
    ScalarShift(TYPE, NAME, shiftSize), origRiskyGrowth(false){
}

/** for reflection */
AdjCreditSpreadRhoParallel::AdjCreditSpreadRhoParallel(): 
    ScalarShift(TYPE, NAME), origRiskyGrowth(false){
}

void AdjCreditSpreadRhoParallel::calculate(TweakGroup*  tweakGroup,
                                           CResults*    results){
    try {
        // get list of names to calculate result for. Remove blanks
        OutputNameArrayConstSP names(this->names(tweakGroup, results));

        if (names->empty()) {
            results->storeNotApplicable(this);            
        } else {
        
            double origPrice;

            for (int idx = 0; idx < names->size(); idx++){
                // store what we want to shift
                setMarketDataName((*names)[idx]);

                // We need to compute a base price for pricing with risky equity before tweaking
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

void AdjCreditSpreadRhoParallel::setRiskyGrowth(bool riskyGrowth) {
    origRiskyGrowth = riskyGrowth;
}

bool AdjCreditSpreadRhoParallel::getRiskyGrowth() const {
    return origRiskyGrowth;
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double AdjCreditSpreadRhoParallel::divisor() const{
    static const string method = "AdjCreditSpreadRhoParallel::divisor";
    try{
        double shiftSize;
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
CClassConstSP AdjCreditSpreadRhoParallel::shiftInterface() const{
    return Shift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP AdjCreditSpreadRhoParallel::restorableShiftInterface() const{
    return RestorableShift::TYPE;
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement the
    Shift interface */
bool AdjCreditSpreadRhoParallel::nameMatches(const OutputName&         name,
                               IObjectConstSP          obj){
    // cast obj to AdjCreditSpreadRhoParallel::Shift and then invoke name method
    const Shift& rhoParallelObj =
        dynamic_cast<const Shift&>(*obj);
    return name.equals(rhoParallelObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement the
    Shift interface */
void AdjCreditSpreadRhoParallel::appendName(OutputNameArray&          namesList,
                              IObjectConstSP          obj){
    // cast obj to AdjCreditSpreadRhoParallel::Shift and then invoke name method
    const Shift& rhoParallelObj =
        dynamic_cast<const Shift&>(*obj);
    OutputNameSP outputName(new OutputName(rhoParallelObj.sensName(this)));
    namesList.push_back(outputName);
}

bool AdjCreditSpreadRhoParallel::shift(IObjectSP obj) {
    // cast obj to AdjCreditSpreadRhoParallel::Shift and then invoke shift method
    Shift& rhoParallelObj =
        dynamic_cast<Shift&>(*obj);
    return rhoParallelObj.sensShift(this);
}

void AdjCreditSpreadRhoParallel::restore(IObjectSP obj) {
    // cast obj to AdjCreditSpreadRhoParallel::Shift and then invoke restore method
    RestorableShift& rhoParallelObj =
        dynamic_cast<RestorableShift&>(*obj);
    rhoParallelObj.sensRestore(this);
}

class AdjCreditSpreadRhoParallelHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new AdjCreditSpreadRhoParallel(AdjCreditSpreadRhoParallel::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new AdjCreditSpreadRhoParallel(shiftSize);
        }
    };
        
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(AdjCreditSpreadRhoParallel, clazz);
        SUPERCLASS(ScalarShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultAdjCreditSpreadRhoParallel);
        FIELD_NO_DESC(origRiskyGrowth);
        FIELD_MAKE_TRANSIENT(origRiskyGrowth);
        // register how to build our sensitivity
        SensitivityFactory::addSens(AdjCreditSpreadRhoParallel::NAME, 
                                    new Factory(), 
                                    new AdjCreditSpreadRhoParallel(AdjCreditSpreadRhoParallel::DEFAULT_SHIFT),
                                    AdjCreditSpreadRhoParallel::Shift::TYPE);
    }

    static IObject* defaultAdjCreditSpreadRhoParallel(){
        return new AdjCreditSpreadRhoParallel();
    }
};

CClassConstSP const AdjCreditSpreadRhoParallel::TYPE = CClass::registerClassLoadMethod(
    "AdjCreditSpreadRhoParallel", typeid(AdjCreditSpreadRhoParallel), AdjCreditSpreadRhoParallelHelper::load);

CClassConstSP const AdjCreditSpreadRhoParallel::Shift::TYPE =
CClass::registerInterfaceLoadMethod(
    "AdjCreditSpreadRhoParallel::Shift", typeid(AdjCreditSpreadRhoParallel::Shift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(AdjCreditSpreadRhoParallel::RestorableShift, clazz);
    EXTENDS(AdjCreditSpreadRhoParallel::Shift);
}
    
CClassConstSP const AdjCreditSpreadRhoParallel::RestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "AdjCreditSpreadRhoParallel::RestorableShift", typeid(AdjCreditSpreadRhoParallel::RestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE
