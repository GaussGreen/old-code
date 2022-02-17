//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : AssetAssetVegaPointwise.cpp
//
//   Description : Controls calculation of Vega pointwise
//
//   Author      : André Segger
//
//   Date        : 24 September 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/AssetVegaPointwise.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/E2CModel.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/ExpiryResult.hpp"
#include "edginc/Results.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/Maths.hpp"


DRLIB_BEGIN_NAMESPACE
AssetVegaPointwise::IShift::~IShift(){} // empty
AssetVegaPointwise::IRestorableShift::~IRestorableShift(){} // empty

/** Sens Control for vega pointwise */
const string AssetVegaPointwise::NAME = "ASSET_VEGA_POINTWISE";
const double AssetVegaPointwise::DEFAULT_SHIFT = 0.001;

/** constructor with explicit shift size */
AssetVegaPointwise::AssetVegaPointwise(double     shiftSize):
    VectorShift(TYPE, NAME, shiftSize){
}

/** constructor with explicit shift size and explicit expiry to
    tweak. Also need to specify all the expiries. Building
    AssetVegaPointwise this way is useful as it allows individal points to
    be tweaked. Note that the supplied expiries are overwritten by the
    vector shift calculate method */
AssetVegaPointwise::AssetVegaPointwise(double                     shiftSize,
                                       const ExpiryConstSP&       expiry,
                                       const ExpiryArrayConstSP&  allExpiries):
    VectorShift(TYPE, NAME, shiftSize){
    this->expiry = expiry; //expiry in parent
    this->cachedExpiries = allExpiries; //expiry in parent
}

/** constructor with explicit shift size and name override (override allows
    a VEGA_POINTWISE calculation to be stored under, eg, VEGA_MATRIX) */
AssetVegaPointwise::AssetVegaPointwise(const string& overrideName,
                                       double     shiftSize):
    VectorShift(TYPE, overrideName, shiftSize){
}

/** for reflection */
AssetVegaPointwise::AssetVegaPointwise():
    VectorShift(TYPE, NAME){
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double AssetVegaPointwise::divisor() const{
const string method = "AssetVegaPointwise::divisor";
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
CClassConstSP AssetVegaPointwise::shiftInterface() const{
    return IShift::TYPE;
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP AssetVegaPointwise::restorableShiftInterface() const{
    return IRestorableShift::TYPE;
}

/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool AssetVegaPointwise::nameMatches(const OutputName&         name,
                                     IObjectConstSP          obj){
    // cast obj to AssetVegaPointwise::Shift and then invoke name method
    const IShift& vegaObj =
        dynamic_cast<const IShift&>(*obj);
    return name.equals(vegaObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's Shift interface */
void AssetVegaPointwise::appendName(OutputNameArray&          namesList,
                               IObjectConstSP          obj){
    // cast obj to AssetVegaPointwise::Shift and then invoke name method
    const IShift& vegaObj =
        dynamic_cast<const IShift&>(*obj);
    OutputNameSP outputName(new OutputName(vegaObj.sensName(this)));
    namesList.push_back(outputName);
}

bool AssetVegaPointwise::shift(IObjectSP obj) {
    // cast obj to AssetVegaPointwise::Shift and then invoke shift method
    IShift& vegaObj =
        dynamic_cast<IShift&>(*obj);
    return vegaObj.sensShift(this);
}

void AssetVegaPointwise::restore(IObjectSP obj) {
    // cast obj to AssetVegaPointwise::Shift and then invoke restore method
    IRestorableShift& vegaObj =
        dynamic_cast<IRestorableShift&>(*obj);
    vegaObj.sensRestore(this);
}

/** The supplied object is queried for the expiries array needed
    for doing vega pointwise and this array is returned. The supplied
    object must implement the AssetVegaPointwise.Shift interface */
IObjectConstSP AssetVegaPointwise::qualifier(IObjectConstSP obj){
    // cast obj to AssetVegaPointwise::Shift and then invoke qualifier method
    const IShift& vegaObj =
        dynamic_cast<const IShift&>(*obj);
    return vegaObj.sensExpiries(this);
}

/** implements a one sided vector derivative for each instance of the
    market data which is sensitive to this SensControl */
void AssetVegaPointwise::calculate(TweakGroup*  tweakGroup,
                                   CResults*    results){
    try{
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
            // see if the instrument has a last sens date method
            LastSensDate* lsd =
                dynamic_cast<LastSensDate*>(tweakGroup->getInstrument());
            DateTime      endDate;
            DateTime      valueDate;

            if (lsd) {
                valueDate = tweakGroup->getInstrument()->getValueDate();
                endDate   = lsd->endDate(this);
            }

            // get base price
            CResults dummyResults;
            double origPrice = calcSensPrice(tweakGroup);
            dummyResults.storePrice(origPrice, "DUMMY");

            for (int idx = 0; idx < names->size(); idx++){
                // store the name of what we want to shift
                setMarketDataName((*names)[idx]);
                /* skip over where result has been calculated already */
                if (!results->exists(this)){
                    try {
                        // now get expiries for which we should tweak
                        cachedExpiries = getExpiries(tweakGroup);
                        // create room for storing the results
                        ExpiryResultArraySP tweaks(new ExpiryResultArray(
                            cachedExpiries->size()));
                        // then loop over the expiries
                        bool expired = false;
                        for (int jdx = 0; jdx < cachedExpiries->size(); jdx++){
                            // to do: ask the instrument when to stop tweaking
                            // store the expiry which we want to tweak
                            expiry = (*cachedExpiries)[jdx];
                            // calculate sens (if not expired)
                            double firstDeriv = 0.0;
                            if (!expired) {
                                firstDeriv = calcOneSidedFirstDeriv(tweakGroup,
                                                                    &dummyResults);
                            }
                            // store result in array
                            (*tweaks)[jdx] = ExpiryResult(expiry, firstDeriv);

                            // do we need to tweak anymore ?
                            if (lsd) {
                                expired = expiry->toDate(valueDate).isGreater(endDate);
                            }
                        }
                        // and store it
                        results->storeGreek(tweaks, this);
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
        results->storeGreek(IObjectSP(new Untweakable(e)), getSensOutputName(),
                            OutputNameSP(new OutputName("")));
    }
}



class AssetVegaPointwiseHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new AssetVegaPointwise(AssetVegaPointwise::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new AssetVegaPointwise(shiftSize);
        }
    };
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(AssetVegaPointwise, clazz);
        SUPERCLASS(VectorShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultAssetVegaPointwise);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(AssetVegaPointwise::NAME,
                                    new Factory(),
                                    new AssetVegaPointwise(AssetVegaPointwise::DEFAULT_SHIFT),
                                    AssetVegaPointwise::IShift::TYPE);
    }

    static IObject* defaultAssetVegaPointwise(){
        return new AssetVegaPointwise();
    }
};

CClassConstSP const AssetVegaPointwise::TYPE = CClass::registerClassLoadMethod(
    "AssetVegaPointwise", typeid(AssetVegaPointwise), AssetVegaPointwiseHelper::load);

CClassConstSP const AssetVegaPointwise::IShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "AssetVegaPointwise::IShift", typeid(AssetVegaPointwise::IShift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(AssetVegaPointwise::IRestorableShift, clazz);
    EXTENDS(AssetVegaPointwise::IShift);
}

CClassConstSP const AssetVegaPointwise::IRestorableShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "AssetVegaPointwise::IRestorableShift",
    typeid(AssetVegaPointwise::IRestorableShift),
    restorableShiftLoad);

DRLIB_END_NAMESPACE
