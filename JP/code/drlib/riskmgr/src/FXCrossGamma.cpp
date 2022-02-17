//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FXCrossGamma.cpp
//
//   Description : FXCrossGamma sensitivity
//
//   Author      : Mark A Robson
//
//   Date        : 16 Mar 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/FXCrossGamma.hpp"
#include "edginc/FXDelta.hpp"
#include "edginc/Delta.hpp"
#include "edginc/ITweakableWithRespectTo.hpp"
#include "edginc/Control.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE

const string FXCrossGamma::NAME = "FX_CROSS_GAMMA";
const double FXCrossGamma::DEFAULT_SHIFT = Delta::DEFAULT_SHIFT;

/** constructor with explicit shift size */
FXCrossGamma::FXCrossGamma(double     shiftSize):
    Sensitivity(TYPE), shiftSize(shiftSize){}

/** for reflection */
FXCrossGamma::FXCrossGamma(): 
    Sensitivity(TYPE), shiftSize(DEFAULT_SHIFT){}

/** Is this sensitivity made using a discrete shift (ie a jump) or a
    an approximately continuous one. Returns false */
bool FXCrossGamma::discreteShift() const{
    return false;
}

/** identifies the name used for storing associated results in the output*/
const string& FXCrossGamma::getSensOutputName() const{
    return NAME;
}

// given  a list of name-pairs for fx cross-gamma, create a list
// of names of underlyings that we would need EQUITY delta for to
// do the cross-gamma calculation - by convention 1st name 
// inside OutputName is the equity, 2nd the fx

static OutputNameArraySP deltaNames(OutputNameArrayConstSP names, int idx) {
    static const string method = "FXCrossGamma::deltaNames";
    try {
        int i;

        OutputNameArraySP deltas(new OutputNameArray(0));

        for (i = 0; i < names->size(); i++) {
            deltas->push_back(OutputNameSP(
                                  new OutputName((*names)[i]->idGet(idx))));
        }

        return deltas;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** override delta calculation to calculate cross gamma */
void FXCrossGamma::calculate(TweakGroup*      tweakGroup,
                             CResults*        results){
    // have to do stock delta before we can do cross gamma
    DeltaSP delta1(new Delta(shiftSize,
                             getModel(),
                             getControl()));
    if (!control->sensitivityRequested(delta1)){
        control->removePacketAfterCalc(Delta::NAME);
        control->removePacketAfterCalc(Delta::SECOND_ORDER_NAME);
    }
    // have to do FX delta before we can do cross gamma
    FXDeltaSP delta2(new FXDelta(shiftSize,
                                 getModel(),
                                 getControl()));
    if (!control->sensitivityRequested(delta2)){
        control->removePacketAfterCalc(FXDelta::NAME);
        control->removePacketAfterCalc(FXDelta::SECOND_ORDER_NAME);
    }
    try {
        if (hasOverrideNames()) {
            delta1->storeOverrideNames(deltaNames(toTweak, 0));
        }
        // set up space for the recording of what shifts we're doing
        shifts = ScalarShiftArray(1);
        // then store delta1 there
        ScalarShiftSP shift1(ScalarShiftSP::dynamicCast((IObjectSP)delta1));
        shifts[0] = shift1;

        delta1->calculate(tweakGroup, results);

        if (hasOverrideNames()) {
            delta2->storeOverrideNames(deltaNames(toTweak, 1));
        }

        // store delta2 
        ScalarShiftSP shift2(ScalarShiftSP::dynamicCast((IObjectSP)delta2));
        shifts[0] = shift2;
        delta2->calculate(tweakGroup, results);
        // then call second order cross derivative calculation
        shifts[0] = shift1;
        shifts.push_back(shift2);
        ScalarShift::calculateCrossDerivative(this,
                                              delta1.get(),
                                              Delta::SECOND_ORDER_NAME,
                                              delta2.get(),
                                              FXDelta::SECOND_ORDER_NAME,
                                              tweakGroup, results);
        shifts.resize(0); // reset shifts
    } catch (exception& e){
        shifts.resize(0); // reset shifts
        throw ModelException(&e,  "FXCrossGamma::calculate");
    }
}

/** When calculating cross gamma, several pricings have to be done. These
    either involve a single shift or a double shift. This routine, called
    upon the Sensitivity returned from Control::getSens(), returns the
    shifts which have been made for the current pricing call */
ScalarShiftArray FXCrossGamma::getComponentShifts() const{
    return shifts;
}

class FXCrossGammaHelper{
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new FXCrossGamma(FXCrossGamma::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new FXCrossGamma(shiftSize);
        }
    };
        
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(FXCrossGamma, clazz);
        SUPERCLASS(Sensitivity);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultFXCrossGamma);
        FIELD(shiftSize, "How big to make the tweak");
        // register how to build our sensitivity
        SensitivityFactory::addSens(FXCrossGamma::NAME, 
                                    new Factory(), 
                                    new FXCrossGamma(FXCrossGamma::DEFAULT_SHIFT),
                                    ITweakableWithRespectTo<Spot>::TYPE);
    }

    static IObject* defaultFXCrossGamma(){
        return new FXCrossGamma();
    }
};

CClassConstSP const FXCrossGamma::TYPE = CClass::registerClassLoadMethod(
    "FXCrossGamma", typeid(FXCrossGamma), FXCrossGammaHelper::load);


DRLIB_END_NAMESPACE
