//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Settlement.cpp
//
//   Description : Abstract settlement interface
//
//   Author      : Andrew J Swain
//
//   Date        : 5 November 2000
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_SETTLEMENT_CPP
#include "edginc/Settlement.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/Addin.hpp"

DRLIB_BEGIN_NAMESPACE

// discount factor between a date & its settlement date
// default concrete implementation
double Settlement::pv(const DateTime& date, YieldCurve *yc) const {
    try {
        return yc->pv(date, settles(date));
    }
    catch (ModelException &e) {
        throw ModelException(&e, "Settlement::pv: Failed");
    }
}

Settlement::Settlement(CClassConstSP clazz): CObject(clazz){}

/** Invoked when Class is 'loaded' */
static void settleLoad(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(Settlement, clazz);
    SUPERCLASS(CObject);
}

CClassConstSP const Settlement::TYPE = CClass::registerClassLoadMethod(
    "Settlement", typeid(Settlement), settleLoad);

Settlement::~Settlement(){}

class SettleAddin: public CObject{
    static CClassConstSP const TYPE;

    SettlementSP    settle;
    DateTimeArraySP tradeDates;     

    /** the 'addin function' */
    static IObjectSP settles(SettleAddin* params){
        static const string routine = "SettleAddin::settles";
        try {
            DateTimeArraySP output(new DateTimeArray(params->tradeDates->size()));
            for (int i = 0; i < params->tradeDates->size(); i++) {
                (*output)[i] = params->settle->settles((*params->tradeDates)[i]);
            }

            return IObjectSP(output);

        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    /** for reflection */
    SettleAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(SettleAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultSettleAddin);
        FIELD(settle, "Settlement");
        FIELD(tradeDates, "Trade Dates");

        Addin::registerClassObjectMethod("SETTLE_DATE",
                                         Addin::RISK,
                                         "Returns settlement date "
                                         "for a trade date",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)settles);
    }

    static IObject* defaultSettleAddin(){
        return new SettleAddin();
    }   
};

CClassConstSP const SettleAddin::TYPE = CClass::registerClassLoadMethod(
    "SettleAddin", typeid(SettleAddin), load);

DRLIB_END_NAMESPACE
