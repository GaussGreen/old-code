//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DeltaProxyNextDay.cpp
//
//   Description : Fund proxy T+1 delta sensitivity
//
//   Author      : Andrew J Swain
//
//   Date        : 8 January 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/DeltaProxyNextDay.hpp"
#include "edginc/DeltaNextDay.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/DeltaProxy.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for DeltaProxyNextDay */
const string DeltaProxyNextDay::NAME = "DELTA_PROXY_NEXT_DAY";
const double DeltaProxyNextDay::DEFAULT_SHIFT = DeltaNextDay::DEFAULT_SHIFT;
const int    DeltaProxyNextDay::DEFAULT_OFFSET = DeltaNextDay::DEFAULT_OFFSET;

/** constructor */
DeltaProxyNextDay::DeltaProxyNextDay(double shift, int offset, HolidaySP hols): 
    T1Delta(TYPE, shift, offset, hols, false) {}

/** for reflection */
DeltaProxyNextDay::DeltaProxyNextDay(): T1Delta(TYPE){}

/** identifies the name used storing associated results in the output */
const string& DeltaProxyNextDay::getSensOutputName() const{
    return NAME;
}

/** From INextDaySensitivity interface */
SensitivitySP DeltaProxyNextDay::requiredSensitivity(
    TweakGroup* tweakGroup) const{
    SensControlPerNameSP delta(new DeltaProxy(deltaShift));
    return Util::requiredSensitivity(this, delta, tweakGroup);
}
    
/** From INextDaySensitivity interface */
void DeltaProxyNextDay::writeResults(const SensitivitySP&          reqdSens,
                                     Results*                      dest,
                                     const UntweakableSP&          untweakable,
                                     const Results*                src) const{
    // can just use simple implementation
    Util::writeResults(this, reqdSens, false, dest, untweakable, src);
}


class DeltaProxyNextDayHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default */
    class Factory: public SensitivityFactory::IDefault {
    public:
        virtual Sensitivity* createDefault(){
            HolidaySP hols(Holiday::weekendsOnly());
            return new DeltaProxyNextDay(DeltaProxyNextDay::DEFAULT_SHIFT,
                                         DeltaProxyNextDay::DEFAULT_OFFSET,
                                         hols);
        }
    };

public:
    /** Invoked when SmallClass is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(DeltaProxyNextDay, clazz);
        SUPERCLASS(T1Delta);
        EMPTY_SHELL_METHOD(defaultDeltaProxyNextDay);
        // register how to build our sensitivity
        SensitivityFactory::addSens(DeltaProxyNextDay::NAME, 
                                    new Factory(), 
                                    new DeltaProxyNextDay(),
                                    0);
    }

    static IObject* defaultDeltaProxyNextDay(){
        return new DeltaProxyNextDay();
    }
};

CClassConstSP const DeltaProxyNextDay::TYPE = CClass::registerClassLoadMethod(
    "DeltaProxyNextDay", typeid(DeltaProxyNextDay), DeltaProxyNextDayHelper::load);

DRLIB_END_NAMESPACE
