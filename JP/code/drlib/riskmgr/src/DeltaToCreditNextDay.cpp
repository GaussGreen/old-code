//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DeltaToCreditNextDay.cpp
//
//   Description : T+1 delta to credit shift
//
//   Author      : André Segger
//
//   Date        : 25 October 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/DeltaToCreditNextDay.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/DeltaToCredit.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for DeltaNextDay */
const string DeltaToCreditNextDay::NAME = "DELTA_TO_CREDIT_NEXT_DAY";
const double DeltaToCreditNextDay::DEFAULT_SHIFT = 0.005;
const int    DeltaToCreditNextDay::DEFAULT_OFFSET = 1;

/** constructor */
DeltaToCreditNextDay::DeltaToCreditNextDay(double shift, 
                                           int offset, HolidaySP hols): 
    T1Delta(TYPE, shift, offset, hols, false) {}

/** for reflection */
DeltaToCreditNextDay::DeltaToCreditNextDay(): T1Delta(TYPE){}

/** identifies the name used storing associated results in the output */
const string& DeltaToCreditNextDay::getSensOutputName() const{
    return NAME;
}

/** From INextDaySensitivity interface */
SensitivitySP DeltaToCreditNextDay::requiredSensitivity(
    TweakGroup* tweakGroup) const{
    SensControlPerNameSP deltaToCredit(new DeltaToCredit(deltaShift));
    return Util::requiredSensitivity(this, deltaToCredit, tweakGroup);
}
   
/** From INextDaySensitivity interface */
void DeltaToCreditNextDay::writeResults(
    const SensitivitySP&          reqdSens,
    Results*                      dest,
    const UntweakableSP&          untweakable,
                                const Results*                src) const{
    Util::writeResults(this, reqdSens, false, dest, untweakable, src);
}

class DeltaToCreditNextDayHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default */
    class Factory: public SensitivityFactory::IDefault {
    public:
        virtual Sensitivity* createDefault(){
            HolidaySP hols(Holiday::weekendsOnly());
            return new DeltaToCreditNextDay(
                DeltaToCreditNextDay::DEFAULT_SHIFT,
                DeltaToCreditNextDay::DEFAULT_OFFSET,
                hols);
        }
    };

public:
    /** Invoked when SmallClass is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(DeltaToCreditNextDay, clazz);
        SUPERCLASS(T1Delta);
        EMPTY_SHELL_METHOD(defaultDeltaToCreditNextDay);
        // register how to build our sensitivity
        SensitivityFactory::addSens(DeltaToCreditNextDay::NAME, 
                                    new Factory(), 
                                    new DeltaToCreditNextDay(),
                                    0);
    }

    static IObject* defaultDeltaToCreditNextDay(){
        return new DeltaToCreditNextDay();
    }
};

CClassConstSP const DeltaToCreditNextDay::TYPE = 
CClass::registerClassLoadMethod(
    "DeltaToCreditNextDay", typeid(DeltaToCreditNextDay), 
    DeltaToCreditNextDayHelper::load);

DRLIB_END_NAMESPACE
