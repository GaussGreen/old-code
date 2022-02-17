//----------------------------------------------------------------------------
//
//   Group       : Convertibles Research
//
//   Filename    : EquityVegaEquivalentParallel.cpp
//
//   Description : ATM equity vega derived from Asset vega 
//
//   Author      : André Segger
//
//   Date        : 04 March 2003
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/EquityVegaEquivalentParallel.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/AssetVegaParallel.hpp"
#include "edginc/Results.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/Control.hpp"


DRLIB_BEGIN_NAMESPACE

/** Sens Control for rho parallel */
const string EquityVegaEquivalentParallel::NAME = "EQUITY_VEGA_EQUIVALENT_PARALLEL";

/** constructor with explicit shift size */
EquityVegaEquivalentParallel::EquityVegaEquivalentParallel(double shiftSize):
    Sensitivity(TYPE), shiftSize(shiftSize){}

/** for reflection */
EquityVegaEquivalentParallel::EquityVegaEquivalentParallel():
    Sensitivity(TYPE){}

/** identifies the name used for storing associated results in the output*/
const string& EquityVegaEquivalentParallel::getSensOutputName() const{
    return NAME;
}

/** Is this sensitivity made using a discrete shift (ie a jump) or a
    an approximately continuous one (return value: false) */
bool EquityVegaEquivalentParallel::discreteShift() const{
    return false;
}

/** implements a one sided scalar derivative for each instance of the
    market data which is sensitive to this SensControl */
void EquityVegaEquivalentParallel::calculate(TweakGroup*  tweakGroup,
                                             CResults*    results) {
    // The implementation here is awful - do not copy from here! */
    try {
        // check that AssetVegaParallel has been requested
        AssetVegaParallelSP assetVega(new AssetVegaParallel(0.0)); // must be 0
        AssetVegaParallelSP assetVegaParallel(
            dynamic_cast<AssetVegaParallel*>(
                control->sensitivityRequested(assetVega).get()));
        
        if (!assetVegaParallel) {
            throw ModelException("EquityVegaEquivalentParallel::calculate",
                    "The calculation of EquityVegaEquivalentParallel requires "
                    "AssetVegaParallel to be requested");
        }
        // calculate asset vega parallel first
        assetVegaParallel->calculateSens(tweakGroup->getModel(),
                                         tweakGroup->getInstrument(),
                                         control,
                                         results);
            
        // get list of names to calculate result for
        OutputNameArrayConstSP names(assetVega->names(tweakGroup));

        if (names->empty()) {
            // store the result
            results->storeNotApplicable(this);
        } else {
            for (int idx = 0; idx < names->size(); idx++) {
                OutputNameConstSP name = (*names)[idx];
                const string& packetName = assetVega->getSensOutputName();
                // retrieve liquidity spread rho results
                if (results->isValidScalarGreek(packetName, name)) {
                    // retrieve liquidity spread rho parallel from results 
                    double assetVegaValue = 
                        results->retrieveScalarGreek(packetName, name);
                    // apply a zero sized perturbation
                    assetVega->findAndShift(IObjectSP::attachToRef(tweakGroup),
                                            name);

                    double equityVega = assetVegaValue * 
                        assetVega->getEquityConversionFactor();
                    OutputNameConstSP outName(
                        new OutputName(assetVega->getEquityName()));
                    // store the result
                    results->storeScalarGreek(equityVega,
                                              getSensOutputName(), outName);
                } else {
                    // this is a complete hack
                    results->storeNotApplicable(this);            
                    return;
                }
            }
        }
    } catch (exception& e){
        results->storeGreek(IObjectSP(new Untweakable(e)), getSensOutputName(),
                            OutputNameSP(new OutputName("")));
    }
}


class EquityVegaEquivalentParallelHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new EquityVegaEquivalentParallel(
                AssetVegaParallel::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new EquityVegaEquivalentParallel(shiftSize);
        }
    };
        
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(EquityVegaEquivalentParallel, clazz);
        SUPERCLASS(Sensitivity);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultEquityVegaEquivalentParallel);
        FIELD(shiftSize, "Not used"); // weak - legacy reasons
        // register how to build our sensitivity
        SensitivityFactory::addSens(EquityVegaEquivalentParallel::NAME, 
                                    new Factory(), 
                                    new EquityVegaEquivalentParallel(
                                        AssetVegaParallel::DEFAULT_SHIFT),
                                    AssetVegaParallel::Shift::TYPE);
    }
    static IObject* defaultEquityVegaEquivalentParallel(){
        return new EquityVegaEquivalentParallel();
    }
};

CClassConstSP const EquityVegaEquivalentParallel::TYPE = 
CClass::registerClassLoadMethod(
    "EquityVegaEquivalentParallel", typeid(EquityVegaEquivalentParallel), 
    EquityVegaEquivalentParallelHelper::load);

DRLIB_END_NAMESPACE
