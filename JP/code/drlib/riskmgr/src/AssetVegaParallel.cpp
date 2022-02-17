
#include "edginc/config.hpp"
#include "edginc/AssetVegaParallel.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/E2CModel.hpp"
#include "edginc/Results.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE
AssetVegaParallel::Shift::~Shift(){} // empty
AssetVegaParallel::RestorableShift::~RestorableShift(){} // empty

/** Sens Control for vega parallel */
const string AssetVegaParallel::NAME = "ASSET_VEGA_PARALLEL";
const double AssetVegaParallel::DEFAULT_SHIFT = 0.001;

/** constructor with explicit shift size */
AssetVegaParallel::AssetVegaParallel(double     shiftSize):
    ScalarShift(TYPE, NAME, shiftSize){
}

/** for reflection */
AssetVegaParallel::AssetVegaParallel(): 
    ScalarShift(TYPE, NAME){
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double AssetVegaParallel::divisor() const{
const string method = "AssetVegaParallel::divisor";
try {
    double shiftSize = getShiftSize();
    if (Maths::isZero(shiftSize)){
       throw ModelException(method, "Shift size is zero");
    }
    return 100.0 * shiftSize;
} catch (exception& e){
    throw ModelException(&e,  method);
    }
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP AssetVegaParallel::shiftInterface() const{
    return Shift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP AssetVegaParallel::restorableShiftInterface() const{
    return RestorableShift::TYPE;
}

/** Returns true if the supplied object matches the supplied name
    for this sensitivity.  The object must implement the
    AssetVegaParallel.Shift interface */
bool AssetVegaParallel::nameMatches(const OutputName&         name,
                                    IObjectConstSP          obj){
    // cast obj to AssetVegaParallel::Shift and then invoke name method
    const Shift& vegaObj =
        dynamic_cast<const Shift&>(*obj);
    return name.equals(vegaObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to
    this sensitivity to the supplied list */
void AssetVegaParallel::appendName(OutputNameArray&          namesList,
                                   IObjectConstSP          obj){
    // cast obj to AssetVegaParallel::Shift and then invoke name method
    const Shift& vegaObj =
        dynamic_cast<const Shift&>(*obj);
    OutputNameSP outputName(new OutputName(vegaObj.sensName(this)));
    namesList.push_back(outputName);
}

bool AssetVegaParallel::shift(IObjectSP obj) {
    // cast obj to AssetVegaParallel::Shift and then invoke shift method
    Shift& vegaObj =
        dynamic_cast<Shift&>(*obj);
    return vegaObj.sensShift(this);
}

void AssetVegaParallel::restore(IObjectSP obj) {
    // cast obj to AssetVegaParallel::Shift and then invoke restore method
    RestorableShift& vegaObj =
        dynamic_cast<RestorableShift&>(*obj);
    vegaObj.sensRestore(this);
}

/** implements a one sided scalar derivative for each instance of the
    market data which is sensitive to this SensControl */
void AssetVegaParallel::calculate(TweakGroup*  tweakGroup,
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
               /* skip over blank names or those where result has
                  been calculated already */
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

void AssetVegaParallel::setEquityConversionFactor(double equityConvFactor)
{
    equityConversionFactor = equityConvFactor;
}

double AssetVegaParallel::getEquityConversionFactor() const
{
    return equityConversionFactor;
}

void AssetVegaParallel::setEquityName(const string& newEquityName)
{
    equityName = newEquityName;
}

string AssetVegaParallel::getEquityName() const
{
    return equityName;
}

class AssetVegaParallelHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new AssetVegaParallel(AssetVegaParallel::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new AssetVegaParallel(shiftSize);
        }
    };
        
public:
    /** Invoked when this class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(AssetVegaParallel, clazz);
        SUPERCLASS(ScalarShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultAssetVegaParallel);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(AssetVegaParallel::NAME, 
                                    new Factory(), 
                                    new AssetVegaParallel(AssetVegaParallel::DEFAULT_SHIFT),
                                    AssetVegaParallel::Shift::TYPE);
    }

    static IObject* defaultAssetVegaParallel(){
        return new AssetVegaParallel();
    }
};

CClassConstSP const AssetVegaParallel::TYPE = CClass::registerClassLoadMethod(
    "AssetVegaParallel", typeid(AssetVegaParallel), AssetVegaParallelHelper::load);

CClassConstSP const AssetVegaParallel::Shift::TYPE =
CClass::registerInterfaceLoadMethod(
    "AssetVegaParallel::Shift", typeid(AssetVegaParallel::Shift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(AssetVegaParallel::RestorableShift, clazz);
    EXTENDS(AssetVegaParallel::Shift);
}
    
CClassConstSP const AssetVegaParallel::RestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "AssetVegaParallel::RestorableShift", typeid(AssetVegaParallel::RestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE
