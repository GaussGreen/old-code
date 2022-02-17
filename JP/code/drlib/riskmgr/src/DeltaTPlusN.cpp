//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DeltaTPlusN.cpp
//
//   Description : T+N delta shift where N is calculated by the instrument and 
//                 is NOT a user input. For instance SPI will calculate delta, 
//                 lag many rebalance dates after today 
//
//   Author      : Ian Stares
//
//   Date        : 23 Nov 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/DeltaTPlusN.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/Delta.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for DeltaTPlusN */
const string DeltaTPlusN::NAME = "DELTA_T_PLUS_N";
const double DeltaTPlusN::DEFAULT_SHIFT = 0.005;
const int    DeltaTPlusN::DEFAULT_OFFSET = 1;

/** constructor */
DeltaTPlusN::DeltaTPlusN(double shift, int offset, HolidaySP hols): 
    T1Delta(TYPE, shift, offset, hols, false) {}

/** for reflection */
DeltaTPlusN::DeltaTPlusN(): T1Delta(TYPE) {}
DeltaTPlusN::DeltaTPlusN(CClassConstSP clazz): T1Delta(clazz) {}

/** identifies the name used storing associated results in the output */
const string& DeltaTPlusN::getSensOutputName() const{
    return NAME;
}

/** From INextDaySensitivity interface */
SensitivitySP DeltaTPlusN::requiredSensitivity(
    TweakGroup* tweakGroup) const{
    SensControlPerNameSP delta(new Delta(deltaShift));
    return Util::requiredSensitivity(this, delta, tweakGroup);
}
   
// override the base class offsetRequired so we can calculate offset
// before then delegating back to base class offsetRequired!
ThetaSP DeltaTPlusN::offsetRequired(const CInstrument* inst) const{
    // use the instrument to calculate how many business days 
    // we should move forward
    const IDeltaTPlusNImnt* viewImnt = dynamic_cast<const IDeltaTPlusNImnt*>(inst);
    if (viewImnt) {
        // do nothing just use whatever offset happens to be there
        offset = viewImnt->deltaOffset(hols.get());
    }

    // now we've set the offset resume normal service 
    // delegate back to the base class for the rest
    return T1Delta::offsetRequired(inst);
}

class DeltaTPlusNHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default */
    class Factory: public SensitivityFactory::IDefault {
    public:
        virtual Sensitivity* createDefault(){
            HolidaySP hols(Holiday::weekendsOnly());
            return new DeltaTPlusN(DeltaTPlusN::DEFAULT_SHIFT,
                                    DeltaTPlusN::DEFAULT_OFFSET,
                                    hols);
        }
    };

public:
    /** Invoked when SmallClass is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(DeltaTPlusN, clazz);
        SUPERCLASS(T1Delta);
        EMPTY_SHELL_METHOD(defaultDeltaTPlusN);
        // register how to build our sensitivity
        SensitivityFactory::addSens(DeltaTPlusN::NAME, 
                                    new Factory(), 
                                    new DeltaTPlusN(),
                                    0);
    }

    static IObject* defaultDeltaTPlusN(){
        return new DeltaTPlusN();
    }
};

CClassConstSP const DeltaTPlusN::IDeltaTPlusNImnt::TYPE =
CClass::registerInterfaceLoadMethod(
    "DeltaTPlusN::IDeltaTPlusNImnt", typeid(DeltaTPlusN::IDeltaTPlusNImnt), 0);

CClassConstSP const DeltaTPlusN::TYPE = CClass::registerClassLoadMethod(
    "DeltaTPlusN", typeid(DeltaTPlusN), DeltaTPlusNHelper::load);

DRLIB_END_NAMESPACE
