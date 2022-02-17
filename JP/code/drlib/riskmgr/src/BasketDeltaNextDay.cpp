//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : BasketDeltaNextDay.cpp
//
//   Description : T+1 delta shift for XCBs
//
//   Author      : Andrew J Swain
//
//   Date        : 20 May 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/T1Delta.hpp"
#include "edginc/BasketDelta.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE

/** XCB T+1 delta shift */
class BasketDeltaNextDay: public T1Delta {
public:
    static CClassConstSP const TYPE;
    const static string NAME;
    const static double DEFAULT_SHIFT;
    const static int    DEFAULT_OFFSET;

    /** identifies the name used storing associated results in the output 
        (Implements pure virtual function in Sensitivity class) */
    const string& getSensOutputName() const {
        return NAME;
    }

    /** From INextDaySensitivity interface */
    virtual SensitivitySP requiredSensitivity(TweakGroup* tweakGroup) const {
        SensControlPerNameSP xcbDelta(new BasketDelta(deltaShift));
        return Util::requiredSensitivity(this, xcbDelta, tweakGroup);
    }
   
    BasketDeltaNextDay(double shift, int offset, HolidaySP hols): 
        T1Delta(TYPE,shift,offset,hols, false) {}

private:
    /** for reflection */
    BasketDeltaNextDay():T1Delta(TYPE) {}

    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default */
    class Factory: public SensitivityFactory::IDefault {
    public:
        virtual Sensitivity* createDefault(){
            HolidaySP hols(Holiday::weekendsOnly());
            return new BasketDeltaNextDay(BasketDeltaNextDay::DEFAULT_SHIFT,
                                          BasketDeltaNextDay::DEFAULT_OFFSET,
                                          hols);
        }
    };

    /** Invoked when SmallClass is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(BasketDeltaNextDay, clazz);
        SUPERCLASS(T1Delta);
        EMPTY_SHELL_METHOD(defaultBasketDeltaNextDay);
        // register how to build our sensitivity
        SensitivityFactory::addSens(BasketDeltaNextDay::NAME, 
                                    new Factory(), 
                                    new BasketDeltaNextDay(),
                                    0);
    }

    static IObject* defaultBasketDeltaNextDay(){
        return new BasketDeltaNextDay();
    }
};

const string BasketDeltaNextDay::NAME = "BASKET_DELTA_NEXT_DAY";
const double BasketDeltaNextDay::DEFAULT_SHIFT = 0.005;
const int    BasketDeltaNextDay::DEFAULT_OFFSET = 1;

CClassConstSP const BasketDeltaNextDay::TYPE = CClass::registerClassLoadMethod(
    "BasketDeltaNextDay", 
    typeid(BasketDeltaNextDay), 
    BasketDeltaNextDay::load);

// to force linker to include file (avoid having header file) */
bool BasketDeltaNextDayLinkIn() {
    return (BasketDeltaNextDay::TYPE != 0);
}

DRLIB_END_NAMESPACE


