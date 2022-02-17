//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SpotPrice.cpp
//
//   Description : SpotPrice sensitivity
//
//   Author      : Stephen Hope
//
//   Date        : 15 May 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SpotPrice.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/ITweakableWithRespectTo.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/BasketDelta.hpp"
#include "edginc/FXDelta.hpp"
#include "edginc/Delta.hpp"
#include "edginc/Results.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/DeltaToCredit.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for SPOT_PRICE */
const string SpotPrice::NAME = "SPOT_PRICE";

/** Is this sensitivity made using a discrete shift (ie a jump) or a
    an approximately continuous one. Returns false */
bool SpotPrice::discreteShift() const{
    return false;
}

/** identifies the name used storing associated results in the output */
const string& SpotPrice::getSensOutputName() const{
    return NAME;
}

static bool isNameRequested(const OutputNameArrayConstSP& requestedNames,
                            const OutputNameSP&           name,
                            CBoolArray&                   spotPriceCalcFlag)
{
    bool isRequested = false;
    if (!requestedNames || requestedNames->size() == 0) {
        isRequested = true;
    } else {
        for(int i=0; i< requestedNames->size();++i) {
            if (name->equals((*requestedNames)[i].get())) {
                isRequested = true;
                spotPriceCalcFlag[i] = true;
                break;
            } 
        }
    }
    return isRequested;
}
                            


/** for reflection */
SpotPrice::SpotPrice(): Sensitivity(TYPE) { }


/** Combines spot prices between results packets (ie does a merge) */
void SpotPrice::addResult(Results*           results,     // (M)
                          const Results*     resultsToAdd,
                          double             scaleFactor) const{
    try{
        results->merge(this->NAME, resultsToAdd);
    } catch (exception& e){
        throw ModelException(e, "SpotPrice::addResult");
    }
}

/** Does the calculation using one specific type of sensitivity
    eg Delta::Shift */
void SpotPrice::subCalculate(const OutputNameArray&        names,
                             IObjectSP                     objToTweak,
                             SensControlPerName&           delta,
                             CBoolArray&                   spotPriceCalcFlag,
                             CResults*                     results){
    for (int idx = 0; idx < names.size(); idx++){
        try {
            if (isNameRequested(toTweak, names[idx],spotPriceCalcFlag)) {
                delta.findAndShift(objToTweak, names[idx]);
                results->storeScalarGreek(delta.getInitialValue(), 
                                          NAME, names[idx]);
            }
        }
        catch (exception& e) {
            results->storeGreek(IObjectSP(new Untweakable(e)),
                                NAME, names[idx]);
        }
    }
}
    
/** Do a DELTA shift only (no pricing) and then store the spot price */
void SpotPrice::calculate(TweakGroup*      tweakGroup,
                          CResults*        results)
{
    static const string method = "SpotPrice::calculate";
    try {
        IObjectSP    objToTweak(IObjectSP::attachToRef(tweakGroup));
        // If you add something here, don't forget to also add it
        // (a) in the "subCalculate" section and
        // (b) in the empty() check for NotApplicable
        Delta         delta(0.0); // must be 0.0 for the shifts
        DeltaToCredit deltaToCredit(0.0);
        BasketDelta   deltabasket(0.0);
        FXDelta       deltaFX(0.0);

        // need to be careful about duplicates etc
        CBoolArray spotPriceCalcFlag;
        if (hasOverrideNames()) {
            spotPriceCalcFlag = BoolArray(toTweak->size());
            for (int i = 0; i < spotPriceCalcFlag.size(); ++i) {
                spotPriceCalcFlag[i] = false;
            }
        }

        OutputNameArrayConstSP names(delta.names(tweakGroup));
        OutputNameArrayConstSP basketNames(deltabasket.names(tweakGroup));
        
        // should possibly strip the names for which we already have a
        // delta output from the e2cNames list
        OutputNameArrayConstSP e2CNames(deltaToCredit.names(tweakGroup));
        OutputNameArrayConstSP fxNames (deltaFX.names(tweakGroup));
        
        if (names->empty()       && e2CNames->empty() &&
            basketNames->empty() && fxNames->empty()) {
            results->storeNotApplicable(this);            
        } else {
            // do spot price for assets implementing Delta::Shift
            subCalculate(*names, objToTweak, delta, 
                         spotPriceCalcFlag, results);
            // do spot price for assets implementing DeltaCredit::Shift
            subCalculate(*e2CNames, objToTweak, deltaToCredit, 
                         spotPriceCalcFlag, results);
            // do spot price for assets implementing DeltaBasket::Shift
            subCalculate(*basketNames, objToTweak, deltabasket, 
                         spotPriceCalcFlag, results);

            // and FX delta for any FX assets
            subCalculate(*fxNames, objToTweak, deltaFX, 
                         spotPriceCalcFlag, results);            
        }

        if (hasOverrideNames()) {
            for (int i = 0;i < spotPriceCalcFlag.size();++i) {
                if (spotPriceCalcFlag[i] == false) {
                    results->storeGreek(
                        IObjectSP(new NotApplicable()), NAME, (*toTweak)[i]);
                }
            }
        }
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

class SpotPriceHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default */
    class Factory: public SensitivityFactory::IDefault {
    public:
        virtual Sensitivity* createDefault(){
            return new SpotPrice();
        }
    };
        
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(SpotPrice, clazz);
        SUPERCLASS(Sensitivity);
        EMPTY_SHELL_METHOD(defaultSpotPrice);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(SpotPrice::NAME, 
                                    new Factory(), 
                                    new SpotPrice(),
                                    ITweakableWithRespectTo<Spot>::TYPE);
    }

    static IObject* defaultSpotPrice(){
        return new SpotPrice();
    }
};


CClassConstSP const SpotPrice::TYPE = CClass::registerClassLoadMethod(
    "SpotPrice", typeid(SpotPrice), SpotPriceHelper::load);


DRLIB_END_NAMESPACE
