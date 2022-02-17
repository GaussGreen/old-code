//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DDeltaDVol.cpp
//
//   Description : DDeltaDVol sensitivity
//
//   Author      : Jay R Blumenstein
//
//   Date        : 15 Apr 2003
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ITweakableWithRespectTo.hpp"
#include "edginc/DDeltaDVol.hpp"
#include "edginc/VegaParallel.hpp"
#include "edginc/VegaParallel2Sided.hpp"
#include "edginc/Delta.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/Maths.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/Control.hpp"
#include "edginc/Results.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/ObjectIteration.hpp"
#include "edginc/IEqVolNamePair.hpp"

DRLIB_BEGIN_NAMESPACE

const string DDeltaDVol::NAME = "D_DELTA_D_VOL";
const double DDeltaDVol::DEFAULT_SHIFT = Delta::DEFAULT_SHIFT;

/** constructor with explicit shift size */
DDeltaDVol::DDeltaDVol(double     shiftSize):
    Sensitivity(TYPE), shiftSize(shiftSize){}

/** for reflection */
DDeltaDVol::DDeltaDVol(): 
    Sensitivity(TYPE), shiftSize(DEFAULT_SHIFT){}

/** Is this sensitivity made using a discrete shift (ie a jump) or a
    an approximately continuous one. Returns false */
bool DDeltaDVol::discreteShift() const{
    return false;
}

/** identifies the name used for storing associated results in the output*/
const string& DDeltaDVol::getSensOutputName() const{
    return NAME;
}

// given  a list of name-pairs for fx cross-gamma, create a list
// of names of underlyings that we would need EQUITY delta for to
// do the cross-gamma calculation - by convention 1st name 
// inside OutputName is the equity, 2nd the fx

static OutputNameArraySP deltaNames(OutputNameArrayConstSP names, int idx) {
    static const string method = "DDeltaDVol::deltaNames";
    try {
        int i;

        OutputNameArraySP deltas(new OutputNameArray(0));

        for (i = 0; i < names->size(); i++) {
            deltas->push_back(OutputNameSP(new OutputName((*names)[i]->idGet(idx))));
        }

        return deltas;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** override delta calculation to calculate cross gamma */
void DDeltaDVol::calculate(TweakGroup*      tweakGroup,
                           CResults*        results) {
    try {
        /** first retrieve equity-vol name pairs */
        MapStringToName namePairs =
            IEqVolNamePair::namePairs(tweakGroup->getInstrumentSP());

        // see if delta is asked for, and if so, use its shift size
        double deltaShift = control->getDeltaShiftSize();
        bool   deltaRequested = !Maths::isZero(deltaShift);
        if (!deltaRequested){
            control->removePacketAfterCalc(Delta::NAME);
            control->removePacketAfterCalc(Delta::SECOND_ORDER_NAME);
        }

        // use supplied Delta shift size, or delta default if not available
        double deltaShiftToUse = !deltaRequested?
            Delta::DEFAULT_SHIFT: deltaShift;
        try {
        
            // do stock delta first
            DeltaSP delta1(new Delta(deltaShiftToUse,
                                     getModel(), getControl()));

            if (hasOverrideNames()) {
                delta1->storeOverrideNames(deltaNames(toTweak, 0));
            }
            // set up space for the recording of what shifts we're doing
            shifts = ScalarShiftArray(1);
            // then store delta1 there
            ScalarShiftSP shift1(ScalarShiftSP::dynamicCast((IObjectSP)delta1));
            shifts[0] = shift1;

            delta1->calculate(tweakGroup, results);

            // see if vega parallel 2 sided is asked for, and if so,
            // use its shift size
            double vega2Shift = control->scalarShiftSize(VegaParallel2Sided::TYPE);
            bool   vega2Requested = !Maths::isZero(vega2Shift);
            if (!vega2Shift){
                control->removePacketAfterCalc(VegaParallel2Sided::NAME);
                control->removePacketAfterCalc(VegaParallel2Sided::SECOND_ORDER_NAME);
            }

            // use supplied vega par 2 sided shift size, or delta
            // default if not available
            double vega2ShiftToUse = !vega2Requested? 
                VegaParallel2Sided::DEFAULT_SHIFT: vega2Shift;

            // have to do vega before we can do DDeltaDVol
            VegaParallel2SidedSP vega(new VegaParallel2Sided(vega2ShiftToUse, 
                                                             getModel(), 
                                                             getControl()));

            if (hasOverrideNames()) {
                vega->storeOverrideNames(deltaNames(toTweak, 1));
            }
            
            // store delta2 
            ScalarShiftSP shift2(ScalarShiftSP::dynamicCast((IObjectSP)vega));
            shifts[0] = shift2;
            vega->calculate(tweakGroup, results);
            // then call second order cross derivative calculation
            shifts[0] = shift1;
            shifts.push_back(shift2);

            // see if we need to bother with anything else
            if (results->isNotApplicable(delta1.get()) ||
                results->isNotApplicable(vega.get())) {
                results->storeNotApplicable(this);
            }
            else {
                ScalarShift::calculateCrossDerivative(
                    this,
                    delta1.get(),
                    Delta::SECOND_ORDER_NAME,
                    vega.get(),
                    VegaParallel2Sided::SECOND_ORDER_NAME,
                    tweakGroup, results,
                    &namePairs);
            }
            shifts.resize(0); // reset shifts
        }
        catch (exception& e){
            OutputNameSP name(new OutputName(""));
            results->storeGreek(IObjectSP(new Untweakable(e)), 
                                getSensOutputName(), name);
        }
    } catch (exception& e){
        throw ModelException(e,  "DDeltaDVol::calculate");
    }
}

/** When calculating cross gamma, several pricings have to be done. These
    either involve a single shift or a double shift. This routine, called
    upon the Sensitivity returned from Control::getSens(), returns the
    shifts which have been made for the current pricing call */
ScalarShiftArray DDeltaDVol::getComponentShifts() const{
    return shifts;
}

class DDeltaDVolHelper{
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new DDeltaDVol(DDeltaDVol::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new DDeltaDVol(shiftSize);
        }
    };
        
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(DDeltaDVol, clazz);
        SUPERCLASS(Sensitivity);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultDDeltaDVol);
        FIELD(shiftSize, "How big to make the tweak");
        // register how to build our sensitivity
        SensitivityFactory::addSens(DDeltaDVol::NAME, 
                                    new Factory(), 
                                    new DDeltaDVol(DDeltaDVol::DEFAULT_SHIFT),
                                    ITweakableWithRespectTo<Spot>::TYPE);
    }

    static IObject* defaultDDeltaDVol(){
        return new DDeltaDVol();
    }
};

CClassConstSP const DDeltaDVol::TYPE = CClass::registerClassLoadMethod(
    "DDeltaDVol", typeid(DDeltaDVol), DDeltaDVolHelper::load);

DRLIB_END_NAMESPACE
