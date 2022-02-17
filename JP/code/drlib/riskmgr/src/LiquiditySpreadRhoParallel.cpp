//----------------------------------------------------------------------------
//
//   Group       : Convertibles Research
//
//   Filename    : LiquiditySpreadRhoParallel.cpp
//
//   Description : liquidity spread rho parallel sensitivity
//
//   Author      : André Segger
//
//   Date        : 04 October 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/LiquiditySpreadRhoParallel.hpp"
#include "edginc/Maths.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/Results.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/E2CModel.hpp"

DRLIB_BEGIN_NAMESPACE
typedef smartPtr<LiquiditySpreadRhoParallel::IRestorableShift> LiquiditySpreadRhoParallelRestorableShiftSP;

// not sure why we have to have these - MSVC gets upset if they're not in the
// header
LiquiditySpreadRhoParallel::Shift::~Shift(){} // empty
LiquiditySpreadRhoParallel::IRestorableShift::~IRestorableShift(){} // empty

const double ONE_BASIS_POINT = 0.0001;

/** Sens Control for rho parallel */
const string LiquiditySpreadRhoParallel::NAME = "LIQUIDITY_SPREAD_RHO_PARALLEL";
const double LiquiditySpreadRhoParallel::DEFAULT_SHIFT = ONE_BASIS_POINT;

/** constructor with explicit shift size */
LiquiditySpreadRhoParallel::LiquiditySpreadRhoParallel(double shiftSize):
    ScalarShift(TYPE, NAME, shiftSize){
}

/** for reflection */
LiquiditySpreadRhoParallel::LiquiditySpreadRhoParallel(): 
    ScalarShift(TYPE, NAME){
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double LiquiditySpreadRhoParallel::divisor() const{
    static const string method = "LiquiditySpreadRhoParallel::divisor";
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
CClassConstSP LiquiditySpreadRhoParallel::shiftInterface() const{
    return Shift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP LiquiditySpreadRhoParallel::restorableShiftInterface() const{
    return IRestorableShift::TYPE;
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement the
    Shift interface */
bool LiquiditySpreadRhoParallel::nameMatches(const OutputName&         name,
                                          IObjectConstSP          obj) {
    // cast obj to RhoParallel::Shift and then invoke name method
    const Shift& rhoParallelObj = 
        dynamic_cast<const Shift&>(*obj);
    return name.equals(rhoParallelObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement the
    Shift interface */
void LiquiditySpreadRhoParallel::appendName(OutputNameArray&          namesList,
                                         IObjectConstSP          obj){
    // cast obj to RhoParallel::Shift and then invoke name method
    const Shift& rhoParallelObj = 
        dynamic_cast<const Shift&>(*obj);
    OutputNameSP outputName(new OutputName(rhoParallelObj.sensName(this)));
    namesList.push_back(outputName);
}

bool LiquiditySpreadRhoParallel::shift(IObjectSP obj) {
    // cast obj to RhoParallel::Shift and then invoke shift method
    Shift& rhoParallelObj = 
        dynamic_cast<Shift&>(*obj);
    return rhoParallelObj.sensShift(this);
}

void LiquiditySpreadRhoParallel::restore(IObjectSP obj) {
    // cast obj to RhoParallel::Shift and then invoke restore method
    IRestorableShift& rhoParallelObj = 
        dynamic_cast<IRestorableShift&>(*obj);
    rhoParallelObj.sensRestore(this);
}

LiquiditySpreadRhoParallel::LiquiditySpreadRhoParallel(double     shiftSize,
                                                       IModel*    model,
                                                       Control*   control):
    ScalarShift(TYPE, NAME, shiftSize) {
    this->algorithm = model;
    this->control   = control;
}

/** implements a one sided scalar derivative for each instance of the
    market data which is sensitive to this SensControl */
void LiquiditySpreadRhoParallel::calculate(TweakGroup*  tweakGroup,
                                           CResults*    results) {
    try {
        // switch to use the correct way of pricing
        if (IE2CModel::TYPE->isInstance(tweakGroup->getModel())) {
            IE2CModel* model = dynamic_cast<IE2CModel*>(tweakGroup->getModel());
            model->setParSpreadPricing(false);
        }

        // get list of names to calculate result for
        OutputNameArrayConstSP names(this->names(tweakGroup, results));

        if (names->empty()) {
            results->storeNotApplicable(this);            
        } else {
            // get base price
            CResults dummyResults;
            double origPrice = calcSensPrice(tweakGroup);
            dummyResults.storePrice(origPrice, "DUMMY");

            for (int idx = 0; idx < names->size(); idx++){
                // store what we want to shift
                setMarketDataName((*names)[idx]);
                /* skip over where result has been calculated already */
                if (!results->exists(this)){
                    try {
                        // calculate sens
                        double firstDeriv = 
                            calcOneSidedFirstDeriv(tweakGroup, &dummyResults);
                        // and store it
                        results->storeScalarGreek(firstDeriv, this);
                    }
                    catch (exception& e) {
                        results->storeGreek(IObjectSP(new Untweakable(e)), this);
                    }
                }
            }
        }

        // switch to use the correct way of pricing
        if (IE2CModel::TYPE->isInstance(tweakGroup->getModel())) {
            IE2CModel* model = dynamic_cast<IE2CModel*>(tweakGroup->getModel());
            model->setParSpreadPricing(true);
        }
    } catch (exception& e){

        // switch back pricing method for the other greeks
        if (IE2CModel::TYPE->isInstance(tweakGroup->getModel())) {
            IE2CModel* model = dynamic_cast<IE2CModel*>(tweakGroup->getModel());
            model->setParSpreadPricing(true);
        }
        throw ModelException(&e,  "ScalarShift::calculate");

    }
}


class LiquiditySpreadRhoParallelHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new LiquiditySpreadRhoParallel(LiquiditySpreadRhoParallel::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new LiquiditySpreadRhoParallel(shiftSize);
        }
    };
        
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(LiquiditySpreadRhoParallel, clazz);
        SUPERCLASS(ScalarShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultRhoParallel);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(LiquiditySpreadRhoParallel::NAME, 
                                    new Factory(), 
                                    new LiquiditySpreadRhoParallel(LiquiditySpreadRhoParallel::DEFAULT_SHIFT),
                                    LiquiditySpreadRhoParallel::Shift::TYPE);
    }

    static IObject* defaultRhoParallel(){
        return new LiquiditySpreadRhoParallel();
    }
};

CClassConstSP const LiquiditySpreadRhoParallel::TYPE = CClass::registerClassLoadMethod(
    "LiquiditySpreadRhoParallel", typeid(LiquiditySpreadRhoParallel), LiquiditySpreadRhoParallelHelper::load);

CClassConstSP const LiquiditySpreadRhoParallel::Shift::TYPE =
CClass::registerInterfaceLoadMethod(
    "LiquiditySpreadRhoParallel::Shift", typeid(LiquiditySpreadRhoParallel::Shift), 0);

static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(LiquiditySpreadRhoParallel::IRestorableShift, clazz);
    EXTENDS(LiquiditySpreadRhoParallel::Shift);
}
    
CClassConstSP const LiquiditySpreadRhoParallel::IRestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "LiquiditySpreadRhoParallel::IRestorableShift", typeid(LiquiditySpreadRhoParallel::IRestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE
