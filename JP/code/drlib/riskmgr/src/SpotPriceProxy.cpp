//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SpotPriceProxy.cpp
//
//   Description : Spot price for fund proxies
//
//   Author      : Andrew J Swain
//
//   Date        : 24 January 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SpotPriceProxy.hpp"
#include "edginc/DeltaProxy.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/Results.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for SPOT_PRICE_PROXY */
const string SpotPriceProxy::NAME = "SPOT_PRICE_PROXY";

/** identifies the name used storing associated results in the output */
const string& SpotPriceProxy::getSensOutputName() const{
    return NAME;
}

/** Is this sensitivity made using a discrete shift (ie a jump) or a
    an approximately continuous one. Returns false */
bool SpotPriceProxy::discreteShift() const{
    return false;
}

/** for reflection */
SpotPriceProxy::SpotPriceProxy(): Sensitivity(TYPE){ }

/** Combines spot prices between results packets (ie does a merge) */
void SpotPriceProxy::addResult(Results*           results,     // (M)
                               const Results*     resultsToAdd,
                               double             scaleFactor) const {
    try {
        results->merge(NAME, resultsToAdd);
    } catch (exception& e){
        throw ModelException(e, "SpotPriceProxy::addResult");
    }
}

/** Do a DELTA shift only (no pricing) and then store the spot price */
void SpotPriceProxy::calculate(TweakGroup*      tweakGroup,
                               CResults*        results)
{
    static const string method = "SpotPriceProxy::calculate";
    try {
        DeltaProxy delta(0.0); // must be 0.0 for the shift
        // get list of names needed to be shifted for shift
        if (hasOverrideNames()) {
            delta.storeOverrideNames(toTweak);
        }

        OutputNameArrayConstSP names(delta.names(tweakGroup));

        if (names->empty()) {
            results->storeNotApplicable(this);            
        } else {
            for (int idx = 0; idx < names->size(); idx++){
                const OutputNameConstSP& name = (*names)[idx];
                try {
                    delta.findAndShift(IObjectSP::attachToRef(tweakGroup),
                                       name);
                    results->storeScalarGreek(delta.getInitialValue(), 
                                              NAME, name);
                }
                catch (exception& e) {
                    results->storeGreek(IObjectSP(new Untweakable(e)), 
                                        NAME, name);
                }
            }
        }

    } catch (exception& e) {
        throw ModelException(&e, method);
    }
}

class SpotPriceProxyHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default */
    class Factory: public SensitivityFactory::IDefault {
    public:
        virtual Sensitivity* createDefault(){
            return new SpotPriceProxy();
        }
    };
        
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(SpotPriceProxy, clazz);
        SUPERCLASS(Sensitivity);
        EMPTY_SHELL_METHOD(defaultSpotPriceProxy);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(SpotPriceProxy::NAME, 
                                    new Factory(), 
                                    new SpotPriceProxy(),
                                    DeltaProxy::IShift::TYPE);
    }

    static IObject* defaultSpotPriceProxy(){
        return new SpotPriceProxy();
    }
};


CClassConstSP const SpotPriceProxy::TYPE = CClass::registerClassLoadMethod(
    "SpotPriceProxy", typeid(SpotPriceProxy), SpotPriceProxyHelper::load);


DRLIB_END_NAMESPACE
