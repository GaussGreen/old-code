//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DeltaNextDay.cpp
//
//   Description : T+1 delta shift
//
//   Author      : Andrew J Swain
//
//   Date        : 16 May 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/DeltaNextDay.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/Delta.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for DeltaNextDay */
const string DeltaNextDay::NAME = "DELTA_NEXT_DAY";
const double DeltaNextDay::DEFAULT_SHIFT = 0.005;
const int    DeltaNextDay::DEFAULT_OFFSET = 1;

/** constructor */
DeltaNextDay::DeltaNextDay(double shift, int offset, HolidaySP hols): 
    T1Delta(TYPE, shift, offset, hols, false) {}

/** for reflection */
DeltaNextDay::DeltaNextDay(): T1Delta(TYPE) {}

/** identifies the name used storing associated results in the output */
const string& DeltaNextDay::getSensOutputName() const{
    return NAME;
}

/** From INextDaySensitivity interface */
SensitivitySP DeltaNextDay::requiredSensitivity(
    TweakGroup* tweakGroup) const{
    SensControlPerNameSP delta(new Delta(deltaShift));
    return Util::requiredSensitivity(this, delta, tweakGroup);
}
   
class DeltaNextDayHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default */
    class Factory: public SensitivityFactory::IDefault {
    public:
        virtual Sensitivity* createDefault(){
            HolidaySP hols(Holiday::weekendsOnly());
            return new DeltaNextDay(DeltaNextDay::DEFAULT_SHIFT,
                                    DeltaNextDay::DEFAULT_OFFSET,
                                    hols);
        }
    };

public:
    /** Invoked when SmallClass is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(DeltaNextDay, clazz);
        SUPERCLASS(T1Delta);
        EMPTY_SHELL_METHOD(defaultDeltaNextDay);
        // register how to build our sensitivity
        SensitivityFactory::addSens(DeltaNextDay::NAME, 
                                    new Factory(), 
                                    new DeltaNextDay(),
                                    0);
    }

    static IObject* defaultDeltaNextDay(){
        return new DeltaNextDay();
    }
};

CClassConstSP const DeltaNextDay::TYPE = CClass::registerClassLoadMethod(
    "DeltaNextDay", typeid(DeltaNextDay), DeltaNextDayHelper::load);

DRLIB_END_NAMESPACE
