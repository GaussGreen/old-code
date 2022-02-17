///////////////////////////////////////////////////////////////////////////////
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DeltaTPlusNHack.hpp
//
//   Description : Delta T+n but implemented to conform to current database "shape" 
//                 since otherwise won't get used for ages
//
//   Date        : Mar 2003
//
///////////////////////////////////////////////////////////////////////////////

#include "edginc/config.hpp"
#include "edginc/DeltaTPlusNHack.hpp"

DRLIB_BEGIN_NAMESPACE

/** constructor */
DeltaTPlusNHack::DeltaTPlusNHack(double shift, int offset, HolidaySP hols): 
    DeltaTPlusN(shift, offset, hols) {}

/** for reflection */
DeltaTPlusNHack::DeltaTPlusNHack(): DeltaTPlusN(TYPE) {}

class DeltaTPlusNHackHelper {
public:
    /** Invoked when this class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(DeltaTPlusNHack, clazz);
        SUPERCLASS(DeltaTPlusN);
        EMPTY_SHELL_METHOD(defaultDeltaTPlusNHack);
    }

    static IObject* defaultDeltaTPlusNHack(){
        return new DeltaTPlusNHack();
    }
};

CClassConstSP const DeltaTPlusNHack::TYPE = CClass::registerClassLoadMethod(
    "DeltaTPlusNHack", typeid(DeltaTPlusNHack), DeltaTPlusNHackHelper::load);


class DeltaTPlus1Hack : public DeltaTPlusN {
public:
    static CClassConstSP const TYPE;
    
private:
    DeltaTPlus1Hack():DeltaTPlusN(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(DeltaTPlus1Hack, clazz);
        SUPERCLASS(DeltaTPlusN);
        EMPTY_SHELL_METHOD(defaultDeltaTPlus1Hack);
    }

    static IObject* defaultDeltaTPlus1Hack(){
        return new DeltaTPlus1Hack();
    }
};

CClassConstSP const DeltaTPlus1Hack::TYPE = 
CClass::registerClassLoadMethod("DeltaTPlus1Hack", typeid(DeltaTPlus1Hack), load);


DRLIB_END_NAMESPACE
